@echo off
REM ============================================
REM Build standalone FastAPI server (Windows)
REM Run from project root
REM ============================================

setlocal ENABLEDELAYEDEXPANSION

REM --------------------------------------------------
REM 1. Resolve path to genetic_code_ncbi.csv
REM --------------------------------------------------
for /f "usebackq delims=" %%i in (`
python - <<EOF
import codonbias, os
print(os.path.join(os.path.dirname(codonbias.__file__), "genetic_code_ncbi.csv"))
EOF
`) do set CSV_PATH=%%i

echo Using genetic_code_ncbi.csv at: %CSV_PATH%

REM --------------------------------------------------
REM 2. Run PyInstaller
REM --------------------------------------------------
pyinstaller ^
  --onedir ^
  --name fastapi_server ^
  --add-data "%CSV_PATH%;codonbias" ^
  --add-data "app\modules\configuration.yaml;modules" ^
  app\api_server.py

if errorlevel 1 (
  echo PyInstaller failed
  exit /b 1
)

REM --------------------------------------------------
REM 3. Copy backend executable into Electron backend
REM --------------------------------------------------
set BACKEND_PATH=ui\DCUB\backend

echo Copying backend executable to: %BACKEND_PATH%

REM Create backend directory if it does not exist
if not exist "%BACKEND_PATH%" (
  mkdir "%BACKEND_PATH%"
)

echo Removing previous executable if exists from: %BACKEND_PATH%\fastapi_server
rmdir /s /q "%BACKEND_PATH%\fastapi_server" 2>nul

echo Copying new executable...
xcopy /E /I /Y "dist\fastapi_server" "%BACKEND_PATH%\fastapi_server"

if errorlevel 1 (
  echo Failed copying FastAPI server
  exit /b 1
)

echo FastAPI server build completed successfully
endlocal
