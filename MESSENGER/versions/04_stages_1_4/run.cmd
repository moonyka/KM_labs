@echo off
setlocal
cd /d "%~dp0..\.."
cmake -S versions/04_stages_1_4 -B build/04_stages_1_4 -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b %errorlevel%
cmake --build build/04_stages_1_4 --config Release -- /m:1 /nr:false
if errorlevel 1 exit /b %errorlevel%
build\04_stages_1_4\Release\messenger.exe fit %*
exit /b %errorlevel%
