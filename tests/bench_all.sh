#!/bin/bash
rm -f bench_results_web
export fg_executable=
./bench_program.sh
export fg_executable=autogrid4
./bench_program.sh
export fg_executable=
