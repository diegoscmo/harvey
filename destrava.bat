@echo off

MODE 20,4

ECHO.
ECHO Destravando pastas..
ECHO.

julia -O3 --check-bounds=no "destrava.jl"