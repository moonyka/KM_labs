@echo off
setlocal
cd /d "%~dp0..\.."
cmake -S versions/03_stages_1_3 -B build/03_stages_1_3 -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b %errorlevel%
cmake --build build/03_stages_1_3 --config Release -- /m:1 /nr:false
if errorlevel 1 exit /b %errorlevel%
ctest --test-dir build/03_stages_1_3 -C Release --output-on-failure
if errorlevel 1 exit /b %errorlevel%
build\03_stages_1_3\Release\messenger.exe doppler %*
exit /b %errorlevel%
