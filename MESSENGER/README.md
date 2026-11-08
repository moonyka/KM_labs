# Практика MESSENGER

Каждая накопительная версия хранится в отдельной папке `versions`. Начните с [разбора первой версии для защиты](versions/01_stage_1/DEFENSE_GUIDE.md).

Общие данные и библиотеки находятся в [релизе common-data-v1](https://github.com/moonyka/KM_labs/releases/tag/common-data-v1). Распакуйте `MESSENGER_common_data.zip` в корень репозитория, сохраняя структуру папок. Нужны Visual Studio 2022 (C++), CMake и Windows x64.

Для запуска первой версии выполните `versions/01_stage_1/run.cmd`.

В архиве имеются исходные уведомления NAIF/CSPICE и IAU SOFA; SHA-256 файлов записаны в `MESSENGER/common_files_sha256.json`.
