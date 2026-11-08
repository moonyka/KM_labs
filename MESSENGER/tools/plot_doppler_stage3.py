from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('--folder', default='results/stage_1_3_teacher')
args = parser.parse_args()
folder = Path(args.folder)
if not folder.is_absolute():
    folder = Path(__file__).resolve().parents[1] / folder
data = np.genfromtxt(folder / 'doppler.csv', delimiter=',', names=True)
eligible = data['eligible_for_fit'] == 1
fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')
for station in np.unique(data['station']):
    mask = eligible & (data['station'] == station)
    ax.scatter(data['tai_offset_s'][mask] / 3600, data['residual_hz'][mask],
               s=4, alpha=.6, label=f'DSS-{int(station)}')
ax.axhline(0, color='black', lw=.6)
ax.set(xlabel='TAI hours since JD 2455643.0', ylabel='Observed minus computed (Hz)',
       title='Two-way Doppler residuals: reference SPK, count time in milliseconds')
ax.grid(alpha=.2)
ax.legend()
fig.savefig(folder / 'doppler_residuals.png', dpi=160)
print(folder / 'doppler_residuals.png')
