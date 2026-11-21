# Запуск в VS Code на Windows

## Подготовка

Нужны Visual Studio 2022 с компонентом «Разработка классических приложений на C++», CMake, VS Code и расширение Microsoft C/C++. Используется MSVC x64, поскольку приложенная CSPICE собрана под него.

После загрузки репозитория распакуйте MESSENGER_common_data.zip из Releases в корень проекта. Должны существовать MESSENGER/data/de442t.bsp, MESSENGER/data/msgr_final.bsp и MESSENGER/vendor/cspice/lib/cspice.lib.

На этом компьютере проект расположен в C:\Users\ionea\Desktop\КомпМат.

## Запуск через F5

1. В VS Code откройте файл рабочей области KM_labs.code-workspace.
2. Слева откройте «Запуск и отладка».
3. Выберите «MESSENGER: копия 1», «копия 2», «копия 3» или «копия 4».
4. Нажмите F5.
5. Выберите режим. Для первой проверки используйте selftest.

Перед запуском CMake автоматически соберёт выбранную копию с отладочной информацией. Вывод появится во встроенном терминале.

| Копия | Режимы |
|---|---|
| 1 | selftest, orbit |
| 2 | selftest, orbit, lighttime |
| 3 | selftest, orbit, lighttime, doppler |
| 4 | selftest, orbit, lighttime, doppler, fit |

Для подготовки к первому этапу откройте versions/01_stage_1/src/main.cpp. Поставьте точку останова на начальном состоянии или вызове Integrator::to, выберите первую копию и orbit.

F9 ставит точку останова, F10 выполняет следующую строку, F11 заходит в функцию. В RelWithDebInfo некоторые локальные переменные могут быть оптимизированы.

Ctrl+Shift+B по умолчанию собирает первую копию. Для другой копии используйте «Терминал → Выполнить задачу» и выберите задачу с её номером.

## Запуск командами

Рабочая папка терминала должна быть MESSENGER:

```powershell
cd MESSENGER
.\versions\01_stage_1\run.cmd
```

Остальные run.cmd находятся в соседних папках копий. Их основные команды — lighttime, doppler и fit.

После сборки:

```powershell
.\build\01_stage_1\Release\messenger.exe orbit --hours 20 --step 60
.\build\02_stages_1_2\Release\messenger.exe lighttime --interval 300
.\build\03_stages_1_3\Release\messenger.exe doppler --trajectory reference
.\build\04_stages_1_4\Release\messenger.exe fit --points 180 --iterations 12
```

Проверка всех копий:

```powershell
.\build_versions.cmd
```

Команды самостоятельно выбирают подпапки в results/corrected. Для отдельного эксперимента задайте --out results/my_experiment.

Флаг --trajectory reference выбирает эталонную SPK для проверки светового времени или доплера. По умолчанию используется интегрированная орбита. Для МНК применяется интегрированная орбита.

Не используйте запуск одного файла через Code Runner: проект требует линковки с CSPICE и, начиная со второго этапа, SOFA.

