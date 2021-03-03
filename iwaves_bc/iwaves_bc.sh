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
makescript=make_scenario_iwaves_bc.py

NUMPROCS=$1

if [ -z "$MPIHOME" ] ; then
    EXEC=$SUN
else
    EXEC="$MPIHOME/bin/mpirun -np $NUMPROCS $SUN"
fi

if [ ! -d $datadir ] ; then
    cp -r $maindatadir $datadir
    echo Creatin input files...
    $PYTHONEXEC scripts/$makescript $datadir
    echo Creating grid...
    $EXEC -g -vv --datadir=$datadir
else
    cp $maindatadir/suntans.dat $datadir/.
fi

echo Running suntans...
$EXEC -s -vv --datadir=$datadir 

