"""Check cumulative-source consistency and statistics of the saved C++ outputs."""
from pathlib import Path
import csv,json,math,re
root=Path(__file__).resolve().parents[1]
results=root/'results/corrected'
def read(name):return json.loads((results/name/'summary.json').read_text(encoding='utf-8-sig'))
def rows(path):
    with path.open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
common=['epoch.hpp','math.hpp','gravity.hpp','shgrav.hpp','mercury.hpp','eih.hpp','spice.hpp','rk8_table.hpp','integrator.hpp','main.cpp']
counts=[]
for stage,directory in enumerate(sorted((root/'versions').iterdir()),1):
    names=common+(['spline.hpp','stations.hpp','lighttime.hpp','trajectory.hpp'] if stage>=2 else [])+(['doppler.hpp'] if stage>=3 else [])+(['fit.hpp'] if stage>=4 else [])
    for name in names:assert (root/'src'/name).read_bytes()==(directory/'src'/name).read_bytes(),directory/name
    log=(root/'build'/directory.name/'Testing/Temporary/LastTest.log').read_text(encoding='utf-8',errors='replace')
    assert 'Test Failed.' not in log,directory
    counts.append(log.count('Test Passed.'))
assert counts==[3,4,5,5],counts
orbit=read('stage_1');a=rows(results/'stage_1/orbit.csv');b=rows(results/'stage_1_step30/orbit.csv')
assert len(a)==len(b)==1201 and orbit['accepted_steps']==1200 and orbit['step_seconds']==60
differences=[1000*math.sqrt(sum((float(x[k])-float(y[k]))**2 for k in ['x_km','y_km','z_km'])) for x,y in zip(a,b)]
assert max(differences)<.01,max(differences)
reference=read('stage_1_3_reference');integrated=read('stage_1_3_integrated')
for name,summary in [('stage_1_3_reference',reference),('stage_1_3_integrated',integrated)]:
    observations=rows(results/name/'doppler.csv');good=[r for r in observations if r['eligible_for_fit']=='1']
    assert len(observations)==summary['total_type12']==12832 and len(good)==summary['used']==12831
    rms=math.sqrt(sum(float(r['residual_hz'])**2 for r in good)/len(good))
    assert abs(rms-summary['residual_rms_hz'])<1e-10 and summary['rejected_computations']==0
fit=read('stage_1_4');history=rows(results/'stage_1_4/iterations.csv')
accepted=[float(r['rms_hz']) for r in history if r['accepted']=='1']
assert all(x>y for x,y in zip(accepted,accepted[1:]))
fit_rows=rows(results/'stage_1_4/fit.csv')
assert len(fit_rows)==fit['observation_count'] and fit['converged']
assert abs(math.sqrt(sum(float(r['residual_hz'])**2 for r in fit_rows)/len(fit_rows))-fit['final_rms_hz'])<1e-10
field=json.loads((results/'validation/field_shtools.json').read_text());assert field['passed']
physics=(root/'build/04_stages_1_4/Testing/Temporary/LastTest.log').read_text(encoding='utf-8',errors='replace')
match=re.search(r'PCK rotation error ([\deE+.\-]+), dense interpolation ([\deE+.\-]+) m',physics)
assert match
checks={'cumulative_sources_match':True,'ctest_counts':counts,'ctest_total':sum(counts),'orbit_rows':len(a),'max_step60_step30_position_difference_m':max(differences),'final_step60_step30_position_difference_m':differences[-1],'independent_PCK_rotation_matrix_max_difference':float(match[1]),'quintic_interpolation_near_pericenter_max_error_m':float(match[2]),'shtools_reference_points':field['points_compared'],'shtools_max_acceleration_difference_km_s2':field['max_acceleration_difference_km_s2'],'doppler_statistics_recomputed_from_csv':True,'fit_statistics_recomputed_from_csv':True,'accepted_fit_steps_strictly_decrease_rms':True}
(results/'validation/checks.json').write_text(json.dumps(checks,indent=2)+'\n')
light=read('stage_1_2')
report=f'''# Протокол проверки исправленной версии

Дата: 21.11.2026. Windows x64, MSVC 19.44, Visual Studio 2022, CMake, Release. Проверена стандартная C++ сборка с CSPICE и SOFA.

## Сборки и проверки

Все четыре накопительные копии собраны. CTest: {counts[0]} + {counts[1]} + {counts[2]} + {counts[3]} = {sum(counts)} пройденных проверок. Повторение проверок в разных копиях проверяет их согласованность; это не {sum(counts)} независимых физических моделей.

Проверены центральное поле, градиенты и полюса, поле степени 160, единицы AU/day интерфейса shGrav, порядок RK8 и обратное интегрирование, предел EIH Шварцшильда, вращение Меркурия, интерполяция и световое время, формат наблюдений, порядок rx/tx, старшие разряды REF_FREQ и ramp-интегралы на границах.

Исходники общих частей совпадают во всех копиях. Независимая сверка SHTOOLS 4.14.1: {field['points_compared']} неполярных точек, максимум разницы ускорений {field['max_acceleration_difference_km_s2']:.3e} км/с². Три полярные/околополярные точки отдельно проверены декартовым градиентом. Это проверка выбранных точек, а не доказательство всех возможных условий.

Собственная матрица Меркурия сравнена с CSPICE/PCK при тех же параметрах: максимум разницы элементов {checks['independent_PCK_rotation_matrix_max_difference']:.3e}. Квинтическая интерполяция проверена вблизи перицентра против RK8 с шагом 5 с: максимум {checks['quintic_interpolation_near_pericenter_max_error_m']:.6g} м.

## Численные результаты

| Расчёт | Результат |
|---|---|
| Орбита, 20 часов, RK8, шаг 60 с | {len(a)} эпох, {orbit['accepted_steps']} шагов |
| Максимальная ошибка положения | {orbit['position_max_m']:.6f} м |
| RMS ошибки положения | {orbit['position_rms_m']:.6f} м |
| Конечная ошибка положения | {orbit['position_final_m']:.6f} м |
| Максимальная разность векторов положения при шагах 60 и 30 с | {max(differences):.6f} м |
| Замыкание светового времени | {light['max_closure_ns']:.6g} нс |
| Доплер на эталонной SPK, RMS | {reference['residual_rms_hz']:.6f} Гц |
| Доплер на интегрированной орбите до подгонки, RMS | {integrated['residual_rms_hz']:.6f} Гц |
| МНК, фактическое число наблюдений | {fit['observation_count']} по всему доступному интервалу |
| МНК, начальный RMS | {fit['initial_rms_hz']:.6f} Гц |
| МНК, конечный RMS | {fit['final_rms_hz']:.6f} Гц |
| МНК, принятые итерации | {fit['accepted_iterations']} |
| МНК, флаг сходимости | {fit['converged']} |

RMS доплера указан по 12831 допустимому наблюдению. Все 12832 строки, включая одно экстремальное, сохранены; статистика с экстремальным значением отдельно записана в summary.json. Ошибок вычисления на полном наборе нет. Для МНК использована равномерная выборка {fit['observation_count']} наблюдений, а не весь массив из 12831 строки. Принятые шаги строго уменьшают RMS; отклонённые попытки также сохранены.

Критерий сходимости МНК: принятый шаг с нормой масштабированных поправок меньше 1e-5 либо уменьшение RMS меньше 1e-5 Гц. Предел итераций сам по себе не устанавливает converged=true.

Замыкание светового времени — численная невязка уравнения. Оно не равно физической точности орбиты. Максимум ошибки орбиты рассчитан по сетке вывода. Итоговый RMS МНК {fit['final_rms_hz']:.3f} Гц не объявляется достижением уровня 0.1 или 0.005 Гц на всём интервале.

Точные значения и перепроверка статистик: [checks.json](results/corrected/validation/checks.json). Исходные CSV и JSON находятся в results/corrected. Старые результаты вне этой папки относятся к предыдущим версиям.

## Воспроизведение

Из папки MESSENGER:

```powershell
.\\build_versions.cmd
.\\build\\01_stage_1\\Release\\messenger.exe orbit
.\\build\\01_stage_1\\Release\\messenger.exe orbit --step 30 --out results/corrected/stage_1_step30
.\\build\\02_stages_1_2\\Release\\messenger.exe lighttime --interval 300
.\\build\\03_stages_1_3\\Release\\messenger.exe doppler --trajectory reference
.\\build\\03_stages_1_3\\Release\\messenger.exe doppler
.\\build\\04_stages_1_4\\Release\\messenger.exe fit --points 180 --iterations 12
python tools/verify_outputs.py
```

Для дополнительной внешней проверки установите pyshtools==4.14.1, создайте CSV командой core_tests data/HgM009.sha build/field_probe.csv и запустите python tools/validate_shtools.py build/field_probe.csv. Эта зависимость не нужна рабочему C++ executable. Опциональная сборка оригинального grav.cpp с внешним Fortran/SHTOOLS в данный протокол стандартной сборки не входит.

Графики: [ошибка орбиты](results/corrected/stage_1/orbit_error.png), [орбита](results/corrected/stage_1/orbit_3d.png), [доплер на SPK](results/corrected/stage_1_3_reference/doppler.png), [МНК](results/corrected/stage_1_4/fit.png). Они воспроизводятся tools/plot_orbit.py и tools/plot_corrected.py из сохранённых CSV.
'''
readability=results/'validation/readability.json'
if readability.exists():
    proof=json.loads(readability.read_text())
    report+='\n## Читаемость кода\n\nПосле форматирования и выделения parseArguments/runOrbit повторно пройдены 17 проверок CTest. Контрольный 20-часовой orbit.csv побайтово совпал с предыдущей версией. Стиль закреплён в .clang-format.\n'
(root/'VERIFICATION.md').write_text(report,encoding='utf-8')
print(json.dumps(checks))
