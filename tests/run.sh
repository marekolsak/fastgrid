#!/bin/bash
if [ -z "$fg_input_file" ]; then
    fg_input_file=../colagen2.gpf
fi
if [ -z "$fg_log_file" ]; then
    fg_log_file=log.glg
fi
if [ -z "$fg_reference_dir" ]; then
    fg_reference_dir=original_win32
fi
if [ -z "$fg_new_dir" ]; then
    fg_new_dir=current
fi
if [ -z "$fg_diff_file" ]; then
    fg_diff_file=diff.txt
fi
if [ -z "$fg_results_file" ]; then
    fg_results_file=results.txt
fi
if [ -z "$fg_test_enabled" ]; then
    fg_test_enabled=n
fi
if [ -z "$fg_executable" ]; then
    fg_executable=fastgrid4
    fg_params="--benchmark --v4"
fi

cd $fg_new_dir

fg_cmd="$fg_executable -p $fg_input_file -l $fg_log_file $fg_params $1 $2 $3 $4 $5 $6 $7 $8 $9"
echo "EXEC: $fg_cmd"
../$fg_cmd
result=$?

if [ $result -eq 0 ]; then
    tail -n 1 $fg_log_file
    cd ..

    if [ "$fg_test_enabled" = "y" ]; then
        fg_cmd="diff -U 0 $fg_reference_dir $fg_new_dir"
        echo "EXEC: $fg_cmd > $fg_diff_file"
        $fg_cmd > $fg_diff_file

        fg_cmd="./testdiff $fg_diff_file $fg_log_file"
        echo "EXEC: $fg_cmd > $fg_results_file"
        $fg_cmd > $fg_results_file

        cat $fg_results_file | head -n 1
        cat $fg_results_file | tail -n 1
    fi
else
    echo "RUNTIME ERROR! Autogrid terminated in an unusual way."
fi
