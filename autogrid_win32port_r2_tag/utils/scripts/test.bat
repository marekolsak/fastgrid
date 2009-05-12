@echo off
if x%1 == x goto ERROR
if x%2 == x goto ERROR
if x%3 == x goto ERROR
if x%4 == x goto ERROR

:MAIN
cd %2
echo cmd: ../autogrid4.exe" -p ../EXAMPLE7.gpf -l EXAMPLE7.glg
"../autogrid4.exe" -p ../EXAMPLE7.gpf -l EXAMPLE7.glg
cd ..
echo cmd: diff -U 0 %1 %2 ^> %3
diff -U 0 %1 %2 > %3
echo cmd: testdiff %3 ^> %4
testdiff %3 > %4
type %4
goto QUIT

:ERROR
echo usage: %0 [reference dir] [autogrid dir] [diff output] [test output]

:QUIT
