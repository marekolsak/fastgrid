#!/bin/bash
if [ -z "$fg_executable" ]; then
    fg_executable=fastgrid4
fi

echo -n "$fg_executable $1 $2: " >> bench_results_web
./rungridmol.sh $1 $1 $1 $2 2>&1|tail -n 1 >> bench_results_web
tail -n 1 bench_results_web