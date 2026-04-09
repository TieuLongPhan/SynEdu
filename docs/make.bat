@ECHO OFF

pushd %~dp0

REM Command file for Jupyter Book documentation

if "%1" == "" goto usage

jupyter-book >NUL 2>NUL
if errorlevel 9009 (
	echo.
	echo.The 'jupyter-book' command was not found.
	echo.Install docs dependencies first:
	echo.  pip install -e ".[docs]"
	echo.
	exit /b 1
)

if "%1" == "html" (
	python ..\script\prepare_doc_notebooks.py
	jupyter-book build .
	goto end
)

if "%1" == "clean" (
	rmdir /s /q _build
	rmdir /s /q talktorials\_generated
	goto end
)

goto usage

:usage
echo.
echo.Usage:
echo.  make.bat html   Build the Jupyter Book
echo.  make.bat clean  Remove build artifacts
goto end

:end
popd
