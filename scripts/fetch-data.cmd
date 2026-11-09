@echo off
rem Downloads the large data sets into data\download (see cmake\FetchData.cmake).
rem   scripts\fetch-data.cmd            fetch what is missing, verify fixed files
rem   scripts\fetch-data.cmd refresh    also update the living series (EOP, Earth PCK)
setlocal
for /f "tokens=2 delims=:." %%c in ('chcp') do set "_cp=%%c"
chcp 65001 >nul
set "_refresh="
if /i "%~1"=="refresh" set "_refresh=-DREFRESH=ON"
cmake %_refresh% -P "%~dp0..\cmake\FetchData.cmake"
set "_rc=%ERRORLEVEL%"
chcp %_cp% >nul
exit /b %_rc%
