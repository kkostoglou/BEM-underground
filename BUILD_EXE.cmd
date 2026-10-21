@echo off
setlocal
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
    echo Create .venv as described in INSTALLATION.md first.
    exit /b 1
)
".venv\Scripts\python.exe" -m pip install -e . pyinstaller==6.22.3
if errorlevel 1 exit /b 1
".venv\Scripts\python.exe" -m PyInstaller --clean --noconfirm BEMUnderground.spec
if errorlevel 1 exit /b 1
echo Built dist\BEMUnderground.exe
