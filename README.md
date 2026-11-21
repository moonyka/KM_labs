# Компьютерная математика: MESSENGER

Четыре накопительные копии: [этап 1](MESSENGER/versions/01_stage_1), [этапы 1–2](MESSENGER/versions/02_stages_1_2), [этапы 1–3](MESSENGER/versions/03_stages_1_3), [этапы 1–4](MESSENGER/versions/04_stages_1_4).

Исправлены динамика EIH, RK8 с шагом 60 с, ориентация Меркурия, работа со временем и согласованность всех копий. Основание: письменные условия и их уточнения до 23.10.2026.

- [Описание проекта и результатов](MESSENGER/README.md)
- [Подготовка к защите первого этапа](MESSENGER/versions/01_stage_1/DEFENSE_GUIDE.md)
- [Сверка с условиями](MESSENGER/CLASS_NOTES_CHECK.md)
- [Происхождение формул и коэффициентов](MESSENGER/MODEL_SOURCES.md)
- [Запуск в VS Code](MESSENGER/RUN_IN_VSCODE.md)
- [Протокол проверок](MESSENGER/VERIFICATION.md)

Откройте KM_labs.code-workspace в VS Code и выберите конфигурацию нужной копии через F5.

Данные и CSPICE/SOFA: скачайте [MESSENGER_common_data.zip](https://github.com/moonyka/KM_labs/releases/download/common-data-v1/MESSENGER_common_data.zip) и распакуйте в корень репозитория. Должны появиться MESSENGER/data и MESSENGER/vendor. Для обычной сборки нужны Visual Studio 2022 с C++ и CMake.

