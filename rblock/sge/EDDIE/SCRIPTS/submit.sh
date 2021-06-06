#!/bin/sh

echo "Running on Eddie..."

WORKINGDIR=/exports/home/cmclean5/STUDIES/HumanSynaptosome/SBM/EDDIE

EXE=$WORKINGDIR/SCRIPTS/execute.sh

nCORES=10

START=1
END=2

name="SBM_20"

chmod +x $EXE

for i in `seq $START $END`
    do
     ###qsub -N $name -l h_rt=00:30:00 -pe sharedmem $nCORES -v nCORES=$nCORES $EXE
     #qsub -N $name -l h_rt=12:00:00 -pe sharedmem $nCORES -v nCORES=$nCORES $EXE
     qsub -N $name -l h_rt=12:00:00 -pe sharedmem $nCORES -l h_vmem=2G -v nCORES=$nCORES $EXE
  done

echo "$0 done!"

exit 0
    
