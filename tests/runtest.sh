#!/bin/bash
cp colagen2.gpf.orig colagen2.gpf
cp macromol_1t60_half_prot.pdbqt.orig macromol_1t60_half_prot.pdbqt

export fg_test_enabled=y
./run.sh $1 $2 $3 $4 $5 $6 $7 $8 $9
export fg_test_enabled=
