"""Create reproducible scientific plots from the C++ calculation output."""
from pathlib import Path
import argparse
import csv
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('folder', nargs='?', default='results/stage_1')
args = parser.parse_args()
folder = Path(args.folder)
with (folder / 'orbit.csv').open(newline='') as file:
    rows = list(csv.DictReader(file))
values = {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}
t = values['elapsed_s'] / 3600
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, layout='constrained')
axes[0].plot(t, values['position_error_m'], color='#155e75', lw=1.6)
axes[0].set_ylabel('Ошибка положения, м')
axes[0].set_title('MESSENGER: сравнение модели с msgr_final.bsp')
axes[1].plot(t, values['velocity_error_m_s'], color='#b45309', lw=1.6)
axes[1].set_ylabel('Ошибка скорости, м/с')
axes[1].set_xlabel('Время от JD 2455646.83 TDB, ч')
for ax in axes:
    ax.grid(alpha=.2)
    ax.set_xlim(0, t[-1])
fig.savefig(folder / 'orbit_error.png', dpi=180)
plt.close(fig)
fig = plt.figure(figsize=(8, 7), layout='constrained')
ax = fig.add_subplot(projection='3d')
ax.plot(values['x_km'], values['y_km'], values['z_km'], label='Модель', lw=1.5)
ax.plot(values['reference_x_km'], values['reference_y_km'], values['reference_z_km'], label='Эфемерида', lw=.9, ls='--')
ax.set(xlabel='X, км', ylabel='Y, км', zlabel='Z, км', title='Орбита относительно Меркурия, J2000')
ax.legend()
ax.set_box_aspect([np.ptp(values[key]) for key in ('x_km', 'y_km', 'z_km')])
fig.savefig(folder / 'orbit_3d.png', dpi=180)
print(json.dumps({'plots': [str(folder / 'orbit_error.png'), str(folder / 'orbit_3d.png')]}))
