"""Plot the corrected C++ outputs without changing observations or statistics."""
from pathlib import Path
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]/'results/corrected'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11})
folder=root/'stage_1_2'
with (folder/'lighttime.csv').open(newline='') as f: rows=list(csv.DictReader(f))
fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True,layout='constrained')
t=[float(r['elapsed_s'])/3600 for r in rows]
for key,label in [('up_s','Земля → аппарат'),('down_s','Аппарат → Земля')]: axes[0].plot(t,[float(r[key]) for r in rows],label=label)
axes[0].set(title='Световое время: интегрированная орбита EIH/RK8',ylabel='Время, с');axes[0].legend()
axes[1].plot(t,[(float(r['up_shapiro_us'])+float(r['down_shapiro_us'])) for r in rows])
axes[1].set(xlabel='Время от начальной эпохи, ч',ylabel='Сумма задержек Шапиро, мкс')
for ax in axes: ax.grid(alpha=.2)
fig.savefig(folder/'lighttime.png',dpi=160);plt.close(fig)
for name,label in [('stage_1_3_reference','эталонная SPK'),('stage_1_3_integrated','интегрированная EIH/RK8')]:
    folder=root/name
    with (folder/'doppler.csv').open(newline='') as f: all_rows=list(csv.DictReader(f))
    rows=[r for r in all_rows if r['eligible_for_fit']=='1'];start=min(float(r['tai_offset_s']) for r in rows)
    stats=json.loads((folder/'summary.json').read_text())
    fig,axes=plt.subplots(2,1,figsize=(11,8),sharex=True,layout='constrained')
    for station in sorted({r['station'] for r in rows},key=int):
        current=[r for r in rows if r['station']==station];t=[(float(r['tai_offset_s'])-start)/3600 for r in current]
        axes[0].scatter(t,[float(r['observed_hz'])/1000 for r in current],s=2,label='DSS-'+station)
        axes[1].scatter(t,[float(r['residual_hz']) for r in current],s=2)
    axes[0].set(title=f'MESSENGER: {label}; RMS = {stats["residual_rms_hz"]:.3f} Гц',ylabel='Наблюдаемый доплер, кГц');axes[0].legend(ncols=3,markerscale=3)
    axes[1].set(xlabel='Время от первого наблюдения, ч',ylabel='Наблюдение − модель, Гц')
    for ax in axes: ax.grid(alpha=.2)
    fig.savefig(folder/'doppler.png',dpi=160);plt.close(fig)
folder=root/'stage_1_4'
with (folder/'fit.csv').open(newline='') as f: rows=list(csv.DictReader(f))
stats=json.loads((folder/'summary.json').read_text());start=min(float(r['tai_offset_s']) for r in rows)
fig,ax=plt.subplots(figsize=(10,5),layout='constrained')
ax.scatter([(float(r['tai_offset_s'])-start)/3600 for r in rows],[float(r['residual_hz']) for r in rows],s=10)
ax.set(title=f'МНК: {len(rows)} наблюдений, RMS = {stats["final_rms_hz"]:.3f} Гц',xlabel='Время от первого наблюдения, ч',ylabel='Невязка после подгонки, Гц');ax.grid(alpha=.2)
fig.savefig(folder/'fit.png',dpi=160);plt.close(fig)
print('Created corrected light-time, Doppler and least-squares plots')
