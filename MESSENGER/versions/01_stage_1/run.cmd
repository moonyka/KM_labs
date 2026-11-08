@echo off
setlocal
cd /d "%~dp0..\.."
cmake -S versions/01_stage_1 -B build/01_stage_1 -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b %errorlevel%
cmake --build build/01_stage_1 --config Release --target messenger -j 4
if errorlevel 1 exit /b %errorlevel%
build\01_stage_1\Release\messenger.exe orbit %*
exit /b %errorlevel%
