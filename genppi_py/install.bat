@echo off
REM Installation script for GenPPI Python interface on Windows

echo GenPPI Python Interface - Windows Installation
echo ===============================================

REM Check Python version and source
python --version 2>nul
if %ERRORLEVEL% NEQ 0 (
    echo Error: Python is not installed or not in PATH.
    echo Please install Python 3.8+ from https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during installation.
    pause
    exit /b 1
)

REM Check if Python version is 3.8 or higher
for /f "tokens=2" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
for /f "tokens=1,2 delims=." %%a in ("%PYTHON_VERSION%") do (
    set MAJOR=%%a
    set MINOR=%%b
)

if %MAJOR% LSS 3 (
    echo Error: Python 3.8 or higher is required.
    echo Current Python version: %PYTHON_VERSION%
    echo Please install Python 3.8+ from https://www.python.org/downloads/
    pause
    exit /b 1
)

if %MAJOR% EQU 3 if %MINOR% LSS 8 (
    echo Error: Python 3.8 or higher is required.
    echo Current Python version: %PYTHON_VERSION%
    echo Please install Python 3.8+ from https://www.python.org/downloads/
    pause
    exit /b 1
)

echo Python version %PYTHON_VERSION% is compatible.

REM Check if using Microsoft Store Python
for /f "tokens=*" %%i in ('where python') do set PYTHON_PATH=%%i
echo Python found at: %PYTHON_PATH%
echo %PYTHON_PATH% | findstr /i "WindowsApps" >nul
if %ERRORLEVEL% EQU 0 (
    echo.
    echo WARNING: Microsoft Store Python detected!
    echo This version may have compatibility issues with genppi-py.
    echo.
    echo RECOMMENDED SOLUTION:
    echo 1. Uninstall Microsoft Store Python
    echo 2. Download Python from https://www.python.org/downloads/
    echo 3. Check "Add Python to PATH" during installation
    echo 4. Restart command prompt and try again
    echo.
    set /p CONTINUE="Do you want to continue anyway? (y/N): "
    if /i not "%CONTINUE%"=="y" (
        echo Installation cancelled.
        pause
        exit /b 1
    )
)

REM Check if pip is available
where pip >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: pip is not installed. Please install pip first.
    pause
    exit /b 1
)

echo.
echo Installing GenPPI Python interface...
echo Note: If installation fails due to compilation errors, you may need:
echo - Microsoft C++ Build Tools
echo - Download from: https://developer.microsoft.com/en-us/windows/downloads/sdk-archive/
echo.

REM Install the package in development mode
pip install -e .

if %ERRORLEVEL% EQU 0 (
    echo.
    echo GenPPI Python interface installed successfully!
    echo.
    echo Usage:
    echo   genppi --help
    echo   genppi -dir samples
    echo.
    echo Note: Use 'genppi' command, not 'genppi.py'
) else (
    echo.
    echo Installation failed!
    echo You may need to install Microsoft C++ Build Tools.
    echo Download from: https://developer.microsoft.com/en-us/windows/downloads/sdk-archive/
)

pause