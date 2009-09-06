#!/bin/bash

run_it()
{
    echo -n "$fg_executable $1 $2: " >> bench_results_web
    ./rungridmol.sh $1 $1 $1 $2 2>&1|tail -n 1 >> bench_results_web
    tail -n 1 bench_results_web
}

bench_one()
{
    export fg_executable=fastgrid4
    run_it $1 $2
    export fg_executable=autogrid4
    run_it $1 $2
    export fg_executable=
}

bench_molecules()
{
if [ "$1" -ge "360" ]; then
    bench_one $1 1000
    bench_one $1 2500
    bench_one $1 5000
    bench_one $1 7500
    bench_one $1 10000
    bench_one $1 15000
fi
    bench_one $1 20000
    bench_one $1 25000
}

bench_grids()
{
#    bench_molecules 10
#    bench_molecules 20
#    bench_molecules 30
#    bench_molecules 40
#    bench_molecules 50
#    bench_molecules 60
#    bench_molecules 70
#    bench_molecules 80
#    bench_molecules 90
#    bench_molecules 100
#    bench_molecules 120
#    bench_molecules 140
#    bench_molecules 160
#    bench_molecules 180
#    bench_molecules 200
#    bench_molecules 233
#    bench_molecules 266
#    bench_molecules 300
#    bench_molecules 350
    bench_molecules 400
    bench_molecules 500
}

#rm -f bench_results_web
bench_grids
