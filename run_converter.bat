@echo off
setlocal
title KFBio Slide Converter - DigitCells

rem Find Python 3.10+: try the py launcher, then python, then python3.
rem A command only counts if it really runs Python and reports 3.10 or newer
rem (the Microsoft Store "python" placeholder prints no version, so it is skipped).
set "PY="
call :try_python py -3
if not defined PY call :try_python python
if not defined PY call :try_python python3
if not defined PY (
    echo.
    echo Python 3.10 or newer was not found. Install it from python.org and tick 'Add Python to PATH' on the first screen.
    echo.
    pause
    exit /b 1
)

rem Install the add-ons once if any are missing.
%PY% -c "import tifffile, numpy, imagecodecs, PIL, tqdm" >nul 2>&1
if errorlevel 1 (
    echo Installing the add-ons the converter needs. This happens only once...
    rem Same list as requirements.txt, for when the launcher was copied without it.
    if exist "%~dp0requirements.txt" (
        %PY% -m pip install -r "%~dp0requirements.txt"
    ) else (
        %PY% -m pip install "tifffile>=2023.7.10" "imagecodecs>=2023.7.10" "numpy>=1.24.0" "pillow>=9.5.0" "tqdm>=4.65.0"
    )
    %PY% -c "import tifffile, numpy, imagecodecs, PIL, tqdm" >nul 2>&1
    if errorlevel 1 (
        echo.
        echo The add-ons could not be installed. Check your internet connection and try again.
        echo.
        pause
        exit /b 1
    )
    echo.
)

%PY% "%~dp0kfbf_converter.py" %*
set "RC=%ERRORLEVEL%"

rem Always pause so results and errors stay visible (double-click or drag-and-drop)
echo.
pause
exit /b %RC%

:try_python
set "PYVER="
for /f "delims=" %%v in ('%* -c "import sys; print('PYOK' if sys.version_info >= (3, 10) else 'PYOLD')" 2^>nul') do set "PYVER=%%v"
if "%PYVER%"=="PYOK" set "PY=%*"
exit /b 0
