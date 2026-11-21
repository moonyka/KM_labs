@echo off
setlocal
cd /d "%~dp0"
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
if errorlevel 1 exit /b %errorlevel%
cmake --build build --config Release -- /m:1 /nr:false
if errorlevel 1 exit /b %errorlevel%
ctest --test-dir build -C Release --output-on-failure
exit /b %errorlevel%
