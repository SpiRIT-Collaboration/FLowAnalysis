#! /bin/bash
#
# DoRPAnalysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh
source runList.sh

##------>>> EDIT HERE 
export MNSFX=rpsim
#export REDO=0

export MNINVERSION=0             ##   <------@@ input
export MNVERSION=50              ##   <------@@ output 
export MNDBVERSION=$MNVERSION
export MDCUT=0.0                 ##   <------@@ mid-rapidity cut
source SetEnvRPSim.sh
export MXEVT=
##<----

 
function flattening() {
    RUN={$MNRNF} SUFX=$MNSFX VER=$MNINVERSION DBVER=$MNDBVERSION root -q -b DoFlattening.C\(11\)
    RUN={$MNRNF} SUFX=$MNSFX VER=$MNINVERSION DBVER=$MNDBVERSION root -q DoFlattening.C\(10\)
}

RUNNUMBER1=(${MNRNF})

function corr(){    ## only the first run
    RUN=${RUNNUMBER1[0]} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root DoRPAnalysis.C\($MXEVT\)
}

function recorr(){    ## only the first run
    export REDO=1
    RUN=${RUNNUMBER1[0]} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root DoRPAnalysis.C\($MXEVT\)
    export REDO=
}


function allcorr(){ ## Go through from the second to the last
    typeset -i I=1
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C  : $LOG
	RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C  : $LOG
	RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C  : $LOG
	RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG &


	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C  : $LOG
	RUN=${RUN} SUFX=$MNSFX VER=$MNINVERSION OVER=$MNVERSION DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG 
	let I++
    done
}

echo ${RUNNUMBER1[0]}
grep function DoRPAnalysis.sh
env | grep MN

echo "type flattenandcorrection <- all full process"
echo "type correction <- correction only"
echo "type flattening <- flatening "

function flattenandcorrection() {
    export REDO=0
    flattening
    corr
    allcorr
}


function correction() {
    corr
    allcorr
}
