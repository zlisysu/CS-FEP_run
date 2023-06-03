ROOT=`pwd`
TARGET='cdk2_28-26'
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