@echo off
set input_file=../EXAMPLE7.gpf
set log_file=EXAMPLE7.glg

if x%1 == x goto ERROR
if x%2 == x goto ERROR
if x%3 == x goto ERROR
if x%4 == x goto ERROR
goto MAIN

:MAIN
cd %2
echo cmd: ../autogrid4.exe -p %input_file% -l %log_file%
"../autogrid4.exe" -p %input_file% -l %log_file% --benchmark --v4
if NOT ERRORLEVEL 0 goto PROGRAMERROR 
tail -n 1 %log_file%
cd ..
echo cmd: diff -U 0 %1 %2 ^> %3
diff -U 0 %1 %2 > %3
echo cmd: testdiff %3 %log_file% ^> %4
testdiff %3 %log_file% > %4
type %4
goto QUIT

:ERROR
echo usage: %0 [reference dir] [autogrid dir] [diff output] [test output]
goto QUIT

:PROGRAMERROR
echo RUNTIME ERROR! Autogrid terminated in an unusual way.
goto QUIT

:QUIT
set input_file=
set log_file=
