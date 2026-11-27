"""Проверить согласованность копий, статистику CSV и сходимость первого этапа."""
from pathlib import Path
import csv
import json
import math

root = Path(__file__).resolve().parents[1]
results = root / 'results/corrected'
common = ['epoch.hpp', 'math.hpp', 'gravity.hpp', 'shgrav.hpp', 'mercury.hpp',
          'mercury_parameters.hpp', 'eih.hpp', 'spice.hpp', 'rk8_table.hpp',
          'integrator.hpp', 'main.cpp']
counts = []
for stage, directory in enumerate(sorted((root / 'versions').iterdir()), 1):
    names = common + (['spline.hpp', 'stations.hpp', 'lighttime.hpp', 'trajectory.hpp']
                      if stage >= 2 else [])
    for name in names:
        assert (root / 'src' / name).read_bytes() == (directory / 'src' / name).read_bytes(), name
    log = root / 'build' / directory.name / 'Testing/Temporary/LastTest.log'
    if log.exists():
        text = log.read_text(errors='replace')
        assert 'Test Failed.' not in text, directory.name
        counts.append(text.count('Test Passed.'))
if counts:
    assert counts == [3, 4, 5, 5], counts

def rows(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        return list(csv.DictReader(stream))

a = rows(results / 'stage_1/orbit.csv')
b = rows(results / 'stage_1_step30/orbit.csv')
assert len(a) == len(b) == 1201
errors = [float(row['position_error_m']) for row in a]
summary = json.loads((results / 'stage_1/summary.json').read_text())
assert errors[0] == 0 and float(a[-1]['elapsed_s']) == 72000
assert summary['accepted_steps'] == 1200 and summary['step_seconds'] == 60
assert abs(max(errors) - summary['position_max_m']) < 1e-9
assert abs(math.sqrt(sum(e*e for e in errors)/len(errors)) - summary['position_rms_m']) < 1e-9
assert abs(errors[-1] - summary['position_final_m']) < 1e-9
assert 'IAU 2009' in summary['mercury_orientation']
delta = max(1000*math.sqrt(sum((float(x[k])-float(y[k]))**2
                              for k in ('x_km', 'y_km', 'z_km')))
            for x, y in zip(a, b))
assert delta < .01, delta
print(json.dumps({'copies_consistent': True, 'ctest_counts': counts,
                  'orbit_rows': len(a), 'statistics_recomputed': True,
                  'max_step60_step30_difference_m': delta}))
