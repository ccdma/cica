#!/bin/bash
#============ PBS Options ============
#QSUB -q pc
#QSUB -W 1:00
#QSUB -A p=1:t=72:c=72
#QSUB -o result.out
#QSUB -m be
#QSUB -M matsuyama.hiroki.24c@st.kyoto-u.ac.jp

#============ Shell Script ============
cd $QSUB_WORKDIR
set -x

# automatically
# export OMP_NUM_THREADS=$QSUB_THREADS

./batch
