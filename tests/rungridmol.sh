#!/bin/bash
sed "s/npts/npts $1 $2 $3/" colagen2.gpf-template| sed "s/dielectric /dielectric -/" > colagen2.gpf
head -n $4 macromol_1t60_half_prot.pdbqt.orig > macromol_1t60_half_prot.pdbqt

./run.sh $5 $6 $7 $8 $9
