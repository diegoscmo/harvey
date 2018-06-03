@echo off
MODE 20,2
start cmd /c julia -O3 --check-bounds=no "roda_todos.jl"
timeout 6
start cmd /c julia -O3 --check-bounds=no "roda_todos.jl"
timeout 6
start cmd /c julia -O3 --check-bounds=no "roda_todos.jl"
timeout 6
start cmd /c julia -O3 --check-bounds=no "roda_todos.jl"