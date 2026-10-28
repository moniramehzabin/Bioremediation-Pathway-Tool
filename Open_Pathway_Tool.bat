@echo off
setlocal
cd /d "%~dp0"
echo Checking the Pathway Tool environment...
where py >nul 2>nul
if not errorlevel 1 (
    py app.py --check >nul 2>nul
    if not errorlevel 1 (
        py app.py
        goto end
    )
)
echo Native Python is not ready. Checking Ubuntu and the bpt environment...
wsl -d Ubuntu-24.04 --cd "%~dp0." bash -lc "for p in ~/miniforge3 ~/miniconda3 ~/anaconda3; do if [ -f \"$p/etc/profile.d/conda.sh\" ]; then source \"$p/etc/profile.d/conda.sh\"; break; fi; done; conda activate bpt && python app.py --check && python app.py"
if errorlevel 1 (
    echo.
    echo Setup is incomplete. Open APP_GUIDE.md for setup instructions.
    echo This launcher does not disable or change Windows security settings.
)
:end
pause
