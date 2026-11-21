"""Independent degree-160 validation against the official SHTOOLS routine.

Usage: python tools/validate_shtools.py build/field_probe.csv --wheel-dir build/reference_python
Alternatively install pyshtools==4.14.1 and omit --wheel-dir.
"""
import argparse,csv,importlib.util,json,math,os
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser()
parser.add_argument('probe',type=Path)
parser.add_argument('--wheel-dir',type=Path)
parser.add_argument('--data',type=Path,default=Path('data/HgM009.sha'))
parser.add_argument('--out',type=Path,default=Path('results/corrected/validation/field_shtools.json'))
args=parser.parse_args()
if args.wheel_dir:
    base=args.wheel_dir.resolve()
    handle=os.add_dll_directory(str(base/'pyshtools.libs'))
    module_path=next((base/'pyshtools').glob('_SHTOOLS*.pyd'))
    spec=importlib.util.spec_from_file_location('_SHTOOLS',module_path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    reference=module.MakeGravGridPoint
else:
    from pyshtools.gravmag import MakeGravGridPoint as reference
cilm=np.zeros((2,161,161));cilm[0,0,0]=1
for line in args.data.read_text().splitlines()[1:]:
    fields=line.replace(',',' ').split();n,m=map(int,fields[:2]);cilm[0,n,m]=float(fields[2]);cilm[1,n,m]=float(fields[3])
errors=[];polar=0
for row in csv.DictReader(args.probe.open(newline='')):
    x,y,z=(float(row[k]) for k in ('x_km','y_km','z_km'))
    r=math.sqrt(x*x+y*y+z*z);rho=math.hypot(x,y)
    if rho/r<1e-5:
        polar+=1
        continue
    latitude=math.atan2(z,rho);longitude=math.atan2(y,x)
    ar,atheta,aphi=reference(cilm,22031.868551,2440.,r,math.degrees(latitude),math.degrees(longitude),lmax=160)
    c,s=math.cos(latitude),math.sin(latitude);cl,sl=math.cos(longitude),math.sin(longitude)
    expected=np.array([ar*c*cl+atheta*s*cl-aphi*sl,ar*c*sl+atheta*s*sl+aphi*cl,ar*s-atheta*c])
    actual=np.array([float(row[k]) for k in ('ax_km_s2','ay_km_s2','az_km_s2')])
    errors.append(float(np.linalg.norm(actual-expected)))
assert len(errors)>=6 and max(errors)<2e-15,(len(errors),max(errors))
result={'reference':'SHTOOLS 4.14.1 MakeGravGridPoint','gravity_degree':160,'GM_km3_s2':22031.868551,'radius_km':2440,'points_compared':len(errors),'polar_points_covered_by_Cartesian_tests':polar,'max_acceleration_difference_km_s2':max(errors),'threshold_km_s2':2e-15,'passed':True}
args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
