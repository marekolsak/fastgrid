#!/bin/bash
echo "Building fastgrid4 ..."
cd ../fastgrid_master/fastgrid/build/
cmake ..
make -j $1
cp fastgrid4 ../../../tests/
cd ../../../tests
