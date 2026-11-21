@echo off
setlocal
call "%~dp0versions\01_stage_1\run.cmd" %*
exit /b %errorlevel%
