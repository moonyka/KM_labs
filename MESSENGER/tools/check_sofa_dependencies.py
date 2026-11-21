from pathlib import Path
import re
import json

folder = Path(__file__).resolve().parents[1] / 'vendor' / 'sofa'
todo = ['xys06a', 'c2ixys', 'era00', 'sp00', 'pom00', 'c2tcio', 'taitt', 'taiutc']
seen = set()
missing = set()
while todo:
    name = todo.pop()
    if name in seen:
        continue
    seen.add(name)
    source = folder / (name + '.c')
    if not source.exists():
        missing.add(name)
        continue
    todo += [m.lower() for m in re.findall(r'\biau([A-Z][A-Za-z0-9]*)\s*\(', source.read_text()) if m.lower() not in seen]
print(json.dumps({'needed': len(seen), 'missing': sorted(missing)}))
