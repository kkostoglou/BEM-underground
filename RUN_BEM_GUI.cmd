@echo off
setlocal
cd /d "%~dp0"
set "PYTHONPATH=%CD%\src"
if exist "%CD%\.conda-env\python.exe" (
    set "PATH=%CD%\.conda-env;%CD%\.conda-env\Library\bin;%CD%\.conda-env\Scripts;%PATH%"
    set "BEU_PYTHON=%CD%\.conda-env\python.exe"
) else if exist "%CD%\.venv\Scripts\python.exe" (
    set "BEU_PYTHON=%CD%\.venv\Scripts\python.exe"
) else (
    echo No project Python environment was found.
    echo Create .conda-env or .venv as described in INSTALLATION.md.
    pause
    exit /b 1
)
"%BEU_PYTHON%" -m python.gui
if errorlevel 1 pause
