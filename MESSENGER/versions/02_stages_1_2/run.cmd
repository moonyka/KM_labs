@echo off
setlocal
cd /d "%~dp0..\.."
cmake -S versions/02_stages_1_2 -B build/02_stages_1_2 -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b %errorlevel%
cmake --build build/02_stages_1_2 --config Release --target messenger -j 4
if errorlevel 1 exit /b %errorlevel%
build\02_stages_1_2\Release\messenger.exe lighttime %*
exit /b %errorlevel%
