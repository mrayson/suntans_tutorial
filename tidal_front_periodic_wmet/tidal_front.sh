#!/bin/sh
########################################################################
#
# Shell script to run a suntans test case.
#
########################################################################

#SUNTANSHOME=../../main
SUNTANSHOME=/home/suntans/code/suntans/main
SUN=$SUNTANSHOME/sun
SUNPLOT=$SUNTANSHOME/sunplot
PYTHONEXEC=python

. $SUNTANSHOME/Makefile.in

maindatadir=rundata
datadir=data
griddir=grids/P96x10
#griddir=grids/P192x32
makescript=make_scenario_tidalfront.py

NUMPROCS=$1

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
    cp $griddir/*.dat $datadir 
    echo Creating input files...
    $PYTHONEXEC scripts/$makescript $datadir
    echo Creating grid...
    $EXEC -g -vvv --datadir=$datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi

echo Running suntans...
$EXEC -s -vv --datadir=$datadir 

