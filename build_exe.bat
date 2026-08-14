@echo off
setlocal

echo ================================================================
echo   Rainfall Frequency Analysis Studio - EXE Builder
echo ================================================================
echo.
echo This script will:
echo   1. Create a local virtual environment (build_env)
echo   2. Install all required packages
echo   3. Package rainfall_app.py into a single standalone .exe
echo      using PyInstaller
echo.
echo IMPORTANT: Just double-click this file, or run it from an
echo ordinary Command Prompt as:   build_exe.bat
echo (Do NOT run it as "python build_exe.bat" - it is a batch
echo  script, not a Python script.)
echo ================================================================
echo.
pause

REM --- Check Python is available ---
where python >nul 2>nul
if errorlevel 1 (
    echo.
    echo [ERROR] Python was not found on PATH.
    echo Please install Python 3.10+ from https://www.python.org/downloads/
    echo and make sure to tick "Add Python to PATH" during installation.
    echo.
    pause
    exit /b 1
)

echo.
echo [1/5] Creating virtual environment "build_env" ...
python -m venv build_env
if errorlevel 1 goto :error

echo.
echo [2/5] Activating virtual environment ...
call build_env\Scripts\activate.bat
if errorlevel 1 goto :error

echo.
echo [3/5] Upgrading pip and installing dependencies ...
python -m pip install --upgrade pip
pip install -r requirements.txt
if errorlevel 1 goto :error

echo.
echo [4/5] Cleaning previous build artifacts (if any) ...
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist RainfallFrequencyAnalysis.spec del /q RainfallFrequencyAnalysis.spec

echo.
echo [5/5] Building the standalone executable with PyInstaller ...
if exist app_icon.ico (
    pyinstaller --noconfirm --onefile --windowed ^
        --name "RainfallFrequencyAnalysis" ^
        --icon=app_icon.ico ^
        rainfall_app.py
) else (
    echo   (no app_icon.ico found next to this script - building without a custom icon)
    pyinstaller --noconfirm --onefile --windowed ^
        --name "RainfallFrequencyAnalysis" ^
        rainfall_app.py
)
if errorlevel 1 goto :error

echo.
echo ================================================================
echo   BUILD COMPLETE
echo ================================================================
echo   Your application is here:
echo   %cd%\dist\RainfallFrequencyAnalysis.exe
echo.
echo   You can copy that single .exe file anywhere and run it -
echo   no Python installation is required on the target machine.
echo ================================================================
echo.
pause
exit /b 0

:error
echo.
echo ================================================================
echo   BUILD FAILED - see the messages above for details.
echo ================================================================
echo.
pause
exit /b 1
