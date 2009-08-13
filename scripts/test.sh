#!/bin/bash
input_file=../colagen2.gpf
log_file=log.glg

if [ "$#" -ge 4 ]; then
    cd $2
    echo cmd: ../fastgrid4 -p $input_file -l $log_file
    ../fastgrid4 -p $input_file -l $log_file --benchmark --v4 $5 $6 $7 $8 $9
    result=$?
    if [ $result -eq 0 ]; then
        tail -n 1 $log_file
        cd ..
        echo "cmd: diff -U 0 $1 $2 > $3"
        diff -U 0 $1 $2 > $3
        echo "cmd: testdiff $3 $log_file > $4"
        ./testdiff $3 $log_file > $4
        cat $4|tail -n 2
    else
	echo "RUNTIME ERROR! Autogrid terminated in an unusual way."
    fi
else
    echo "usage: $0 [reference dir] [autogrid dir] [diff output] [test output]"
fi
