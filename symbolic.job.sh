#!/bin/bash
#============ PBS Options ============
#QSUB -q pc
#QSUB -W 24:00
#QSUB -A p=1:t=72:c=72

#============ Shell Script ============
cd $QSUB_WORKDIR
set -x

# automatically
# export OMP_NUM_THREADS=$QSUB_THREADS

./symbolic.out
