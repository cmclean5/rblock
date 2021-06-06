#!/bin/sh
#$ -e /exports/eddie/scratch/cmclean5
#$ -o /exports/eddie/scratch/cmclean5

SUBDIR=$JOB_ID
echo "SUBDIR is $SUBDIR"

EXECDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/SBM
SCRIPTDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/SBM/EDDIE/SCRIPTS
DATADIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/SBM

DATA=datasets.tar.gz
SEED=$JOB_ID
CORES=$nCORES
NITS=10
NTWRK="PPI_PSP"
K=46

cd $TMPDIR
echo "WORKING DIR " $TMPDIR

cp -r $DATADIR/$DATA .
cp -r $SCRIPTDIR/PSP.R .

# initiallise environment module
. /etc/profile.d/modules.sh

#load module R
module load igmm/apps/R/3.6.3 #for rblock library

tar -xzvf $DATA

time Rscript PSP.R $SEED $NITS $CORES $NTWRK $K

OUTDIR=$EXECDIR/EDDIE/RESULTS/$SUBDIR

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

cp -v *.RDS $OUTDIR
cp -v *.csv $OUTDIR

echo "$0 done!"

exit 0
