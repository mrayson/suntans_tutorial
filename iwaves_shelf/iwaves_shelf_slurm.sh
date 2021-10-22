#!/bin/bash --login
#
#SBATCH --account=pawsey0106
#SBATCH --time=01:00:00
#SBATCH --export=NONE
#SBATCH --partition=workq
##SBATCH --partition=debugq
#SBATCH --output=LOGS/Browse2d-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=matt.rayson@uwa.edu.au
##SBATCH --nodes=1
#SBATCH --ntasks=4

########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

#SUNTANSHOME=../../main
SUNTANSHOME=$MYGROUP/testing/suntans/main
SUN=$SUNTANSHOME/sun

## Local python environment
#module use ~/code/modulefiles
#module load anaconda-python/3.6.0
#
#PYTHONEXEC="srun -u --export=all -n 1 -d 1 python"

# Containerised python environment
module load singularity

containerDir=/group/pawsey0106/mrayson/singularity
export containerImage=$containerDir/python_sfoda006.sif
PYTHONEXEC="srun -u --export=all -n $SLURM_NTASKS singularity exec $containerImage python -u"

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
makescript=make_scenario_iwaves.py

# Local output folder
#datadir=data

# Scratch outputfolder
datadir=$MYSCRATCH/SUNTANS_test

NUMPROCS=$SLURM_NTASKS

EXEC="srun -n $NUMPROCS $SUN"

#if [ ! -d $datadir ] ; then
#    cp -r $maindatadir $datadir
#    echo Creatin input files...
#    $PYTHONEXEC scripts/$makescript $datadir
#    echo Creating grid...
#    $EXEC -g -vvv --datadir=$datadir
#else
#    cp $maindatadir/suntans.dat $datadir/.
#fi

cp -r $maindatadir $datadir
echo Creating input files...
$PYTHONEXEC scripts/$makescript $datadir
echo Creating grid...
$EXEC -g -vvv --datadir=$datadir

echo Running suntans...
$EXEC -s -vvv --datadir=$datadir 

