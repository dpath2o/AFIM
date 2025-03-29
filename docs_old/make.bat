@ECHO OFF

REM Minimal make.bat for building Sphinx documentation on Windows

set SPHINXBUILD=sphinx-build
set SOURCEDIR=.
set BUILDDIR=_build

if "%1"=="" (
    %SPHINXBUILD% -M help %SOURCEDIR% %BUILDDIR%
    goto end
)

%SPHINXBUILD% -M %1 %SOURCEDIR% %BUILDDIR%

:end
