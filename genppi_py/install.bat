@echo off
REM Installation script for GenPPI Python interface on Windows

REM Check if pip is available
where pip >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo Error: pip is not installed. Please install pip first.
    exit /b 1
)

REM Install the package in development mode
pip install -e .

echo GenPPI Python interface installed successfully!
echo You can now use the 'genppi' command from the command line.