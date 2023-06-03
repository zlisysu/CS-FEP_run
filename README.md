# System Requirements
## Hardware requirements
CS-FEP requires a standard computer with GPUs like NVIDIA GeForce RTX 2080 Ti or Tesla V100-SXM3-32GB.
## Software requirements
### OS Requirements
This package is supported for *Linux*. The package has been tested on the following systems:
+ Linux: CentOs 7

### Simulation Software
+ Modified Amber with version >= 18. The source code file ti.F90 in AMBERHOME/src/pmemd/src was modified in order to avoid reordering the atoms of ligand.
- To modify Amber software, patch files are provided for amber18 and amber 20. 
    - For amber18: use the `ti.F90.patch` and `patch.sh` in amber18_patch
    - For amber20: use the `ti.F90.patch` and `patch.sh` in amber20_patch
    
The cotent of `patch.sh`
```sh 
#!/bin/bash
# cd the path ($AMBERHOME//src/pmemd/src) that the ti.F90 file exists.
patch -p0 < ti.F90.patch
```    
The installation guide of Amber: http://ambermd.org/Installation.php
- Typical install time of Amber on a "normal" desktop computer: about 40 mins.
### Analysis Code Dependencies
- The analysis code depends on the Python scientific stack. 
```
numpy
pandas
pymbar
matplotlib
```
- Typical install time of python and required packages on a "normal" desktop computer: about 20 mins.


# Running simulation
- Use this script to run the MD simulation
`run.sh`
```sh
#!/bin/bash
ROOT=`pwd`
TARGET='cdk2_28-26_demo'
WINS='0.00 0.05 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 0.95 1.00'
SIDES='complex/cM2A complex/cM2B ligands/lM2A ligands/lM2B'
_CUDA_VISIBLE_DEVICES_=0
for side in $SIDES;do
    cd $ROOT/$TARGET/$side
    for win in $WINS;do
        cd $win
        sh submit.sh $_CUDA_VISIBLE_DEVICES_
    done
done
```
- Expected run time for demo's MD simulation on a "normal" desktop computer(with one gpu like Tesla V100-SXM3-32GB): 6-8 hours


# Analysis
- The demo's simulation out was provided in the cdk2_28-26.tar.bz2, use the following command to extract the data:
```sh
tar -zxvf cdk2_28-26.tar.bz2
```
- Use the in-house free energy calculation python script: `analysis_free_energy.py`
- `aly.sh`
```sh
#!/bin/bash
ROOT=`pwd`
FEP_PATH="$ROOT/cdk2_28-26/FEP"
mkdir -p analysis_out
cd analysis_out
python $ROOT/analysis_free_energy.py -d ${FEP_PATH} -f 0.25 -o fe_out.csv
cd ..
```
- Expected run time for demo's analysis on a "normal" desktop computer: 3-5 mins
