#!/bin/bash --login
#
#SBATCH --account=pawsey0106
#SBATCH --time=00:10:00
#SBATCH --export=NONE
##SBATCH --partition=workq
#SBATCH --partition=debugq
#SBATCH --output=LOGS/Browse2d-%j.out
##SBATCH --mail-type=ALL
##SBATCH --mail-user=matt.rayson@uwa.edu.au
#SBATCH --nodes=1
#SBATCH --ntasks=24

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
export containerImage=$containerDir/python_sfoda007.sif
PYTHONEXEC="srun -u --export=all -n 1 singularity exec $containerImage python -u"

. $SUNTANSHOME/Makefile.in

griddir=grids/P200x250
maindatadir=rundata
makescript=make_scenario_tidalfront.py
plotscript=

# Local output folder
#datadir=data

# Scratch outputfolder
datadir=$MYSCRATCH/SUNTANS_test/iwaves_tidalfront

if [ ! -d $datadir ] ; then
    mkdir $datadir
fi

NUMPROCS=${SLURM_NTASKS}
echo 'Running on NUMPROCS '$NUMPROCS

EXEC="srun -n $NUMPROCS -N ${SLURM_JOB_NUM_NODES} $SUN"

#if [ ! -d $datadir ] ; then
#    cp -r $maindatadir $datadir
#    echo Creatin input files...
#    $PYTHONEXEC scripts/$makescript $datadir
#    echo Creating grid...
#    $EXEC -g -vvv --datadir=$datadir
#else
#    cp $maindatadir/suntans.dat $datadir/.
#fi

# cleanup the working directory
rm $datadir/*.nc 

cp $maindatadir/* $datadir
cp $griddir/*.dat $datadir 
echo Creating input files...
$PYTHONEXEC scripts/$makescript $datadir
echo Creating grid...
$EXEC -g -vvv --datadir=$datadir

echo Running suntans...
$EXEC -s -vvv --datadir=$datadir 

#echo Plotting solution...
#$PYTHONEXEC scripts/plot_vslice_snaps.py "$datadir/IWaveRidge_00*.nc" $datadir/iwave_ridge_gls3
#$PYTHONEXEC scripts/plot_vslice_nuv.py "$datadir/IWaveRidge_00*.nc" $datadir/iwave_ridge_nu_gls3
#
