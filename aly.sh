#!/bin/bash
ROOT=`pwd`
FEP_PATH="$ROOT/cdk2_28-26/FEP"
mkdir -p analysis_out
cd analysis_out
python $ROOT/analysis_free_energy.py -d ${FEP_PATH} -f 0.25 -o fe_out.csv
cd ..
