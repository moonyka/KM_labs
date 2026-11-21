"""Refresh separate cumulative source copies; share input data and libraries."""
from pathlib import Path
import shutil

root = Path(__file__).resolve().parents[1]
source = root / 'src'
cmake = (root / 'CMakeLists.txt').read_text(encoding='utf-8')
cmake = cmake.replace('${CMAKE_CURRENT_SOURCE_DIR}/vendor/cspice', '${CMAKE_CURRENT_SOURCE_DIR}/../../vendor/cspice')
cmake = cmake.replace('"vendor/sofa/*.c"', '"../../vendor/sofa/*.c"').replace('PUBLIC vendor/sofa', 'PUBLIC ../../vendor/sofa')

def remove_fit(text):
    replacements = {
        '#include "fit.hpp"\n': '',
        'command=="fit"?"results/stage_1_4":': '',
        ',fitHours=0': '',
        ',points=180,iterations=3': '',
        'else if(key=="--points")points=std::stoi(value);else if(key=="--iterations")iterations=std::stoi(value);else if(key=="--fit-hours")fitHours=std::stod(value);': '',
        '||command=="fit"': '',
        'orbit, lighttime, doppler, fit, selftest': 'orbit, lighttime, doppler, selftest',
        '&&points>=12&&iterations>=1&&fitHours>=0': '',
        '    if(command=="fit"){runFit(data,out,points,iterations,fitHours);kclear_c();return 0;}\n': '',
    }
    for before, after in replacements.items():
        text = text.replace(before, after)
    return text

def remove_doppler(text):
    replacements = {
        '#include "doppler.hpp"\n': '',
        'command=="doppler"?"results/stage_1_3":': '',
        ',stride=1': '',
        'else if(key=="--stride")stride=std::stoi(value);': '',
        '||command=="doppler"': '',
        'orbit, lighttime, doppler, selftest': 'orbit, lighttime, selftest',
        '&&stride>0': '',
        '    if(command=="doppler"){runDoppler(data,out,stride);kclear_c();return 0;}\n': '',
    }
    for before, after in replacements.items():
        text = text.replace(before, after)
    return text

names = ['01_stage_1', '02_stages_1_2', '03_stages_1_3', '04_stages_1_4']
commands = ['orbit', 'lighttime', 'doppler', 'fit']
for stage, name in enumerate(names, 1):
    folder = root / 'versions' / name
    folder.mkdir(parents=True, exist_ok=True)
    if stage > 1:
        (folder / 'src').mkdir(exist_ok=True)
        for file in source.iterdir():
            if file.name == 'fit.hpp' and stage < 4:
                continue
            if file.name == 'doppler.hpp' and stage < 3:
                continue
            text = file.read_text(encoding='utf-8')
            if file.name == 'main.cpp':
                if stage < 4:
                    text = remove_fit(text)
                if stage < 3:
                    text = remove_doppler(text)
            (folder / 'src' / file.name).write_text(text, encoding='utf-8')
        (folder / 'CMakeLists.txt').write_text(cmake, encoding='utf-8')
    (folder / 'README.md').write_text(
        f'# Накопительная версия: этапы 1–{stage}\n\n'
        f'Запустите `run.cmd`. Команда по умолчанию: `{commands[stage-1]}`. '
        'Данные находятся в общей папке `../../data`, библиотеки — в `../../vendor`. '
        'Основное описание методов и ограничений: `../../README.md`.\n', encoding='utf-8')
    script = f'''@echo off
setlocal
cd /d "%~dp0..\\.."
cmake -S versions/{name} -B build/{name} -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b %errorlevel%
cmake --build build/{name} --config Release --target messenger -j 4
if errorlevel 1 exit /b %errorlevel%
build\\{name}\\Release\\messenger.exe {commands[stage-1]} %*
exit /b %errorlevel%
'''
    (folder / 'run.cmd').write_text(script, encoding='ascii', newline='\r\n')
print('Saved four cumulative source versions')
