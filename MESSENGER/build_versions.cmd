@echo off
setlocal
cd /d "%~dp0"
for %%v in (01_stage_1 02_stages_1_2 03_stages_1_3 04_stages_1_4) do (
  cmake -S versions/%%v -B build/%%v -G "Visual Studio 17 2022" -A x64
  if errorlevel 1 exit /b 1
  cmake --build build/%%v --config Release -- /m:1 /nr:false
  if errorlevel 1 exit /b 1
  ctest --test-dir build/%%v -C Release --output-on-failure
  if errorlevel 1 exit /b 1
)
