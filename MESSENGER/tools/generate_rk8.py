"""Generate DP8 coefficients from the source linked in the practice condition."""
import ast
import shutil
import subprocess
import hashlib
import urllib.request
from fractions import Fraction
from pathlib import Path
root=Path(__file__).resolve().parents[1]
source=root/'build/reference/RK_Butcher_Table_Dictionary.py'
source.parent.mkdir(parents=True,exist_ok=True)
url='https://raw.githubusercontent.com/zachetienne/nrpytutorial/6c83fe6ec1134bf2941e3ce527942d8852920c14/MoLtimestepping/RK_Butcher_Table_Dictionary.py'
expected='4c914c27966e9b1f13d073ffde6695f01c8241c00094e21bfc917e6e2b6f676f'
if not source.exists():
    local=root/'sources/RK_Butcher_Table_Dictionary.py'
    source.write_bytes(local.read_bytes() if local.exists() else urllib.request.urlopen(url).read())
assert hashlib.sha256(source.read_bytes()).hexdigest()==expected,'Upstream coefficient source changed' 
tree=ast.parse(source.read_text(encoding='utf-8'))
assignment=next(n for n in tree.body if isinstance(n,ast.Assign) and isinstance(n.targets[0],ast.Subscript) and isinstance(n.targets[0].slice,ast.Constant) and n.targets[0].slice.value=='DP8')
def number(n):
    if isinstance(n,ast.Constant): return n.value if isinstance(n.value,str) else Fraction(n.value)
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub): return -number(n.operand)
    if isinstance(n,(ast.List,ast.Tuple)): return [number(x) for x in n.elts]
    if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr in ('Rational','sympify'): return Fraction(*(number(x) for x in n.args))
    raise ValueError(ast.dump(n))
table,order=number(assignment.value)
assert order==8 and len(table)==14
def cpp(x): return f'{x.numerator}./{x.denominator}.' if x.denominator!=1 else f'{x.numerator}.'
rows=['{'+','.join(cpp(x) for x in row[1:])+'}' for row in table[:-1]]
text='#pragma once\ninline constexpr int RK_STAGES=13;\n'
text+='inline constexpr double RK_C[13]={'+','.join(cpp(row[0]) for row in table[:-1])+'};\n'
text+='inline constexpr double RK_A[13][13]={\n'+',\n'.join(rows)+'\n};\n'
text+='inline constexpr double RK_B[13]={'+','.join(cpp(x) for x in table[-1][1:])+'};\n'
(root/'src/rk8_table.hpp').write_text(text,encoding='utf-8')
print('Generated 13-stage DP8 tableau')

formatter=shutil.which('clang-format')
if not formatter and Path('C:/Program Files/LLVM/bin/clang-format.exe').exists():
    formatter='C:/Program Files/LLVM/bin/clang-format.exe'
if formatter:
    subprocess.run([formatter,'-i','-style=file',str(root/'src/rk8_table.hpp')],check=True)
