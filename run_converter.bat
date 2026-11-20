@echo off
setlocal
title KFBio Slide Converter - DigitCells

python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python is not installed or not in PATH.
    echo Please install Python 3.10+ and check 'Add Python to PATH'.
    pause
    exit /b 1
)

python "%~dp0kfbf_converter.py" %*
set "RC=%ERRORLEVEL%"

rem Always pause so results and errors stay visible (double-click or drag-and-drop)
echo.
pause
exit /b %RC%
