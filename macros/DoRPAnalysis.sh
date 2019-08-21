#! /bin/bash
#
# DoRPAnalysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh
source runList.sh

##------>>> EDIT HERE 
export MNRNF=$RNF108
export MNDBVERSION=37
export MNVERSION=37.1
export MNDB=BTt
##<----

function flattening() {
    RUN={$MNRNF} VER=$MNDBVERSION SUFX=$MNDB root -q -b DoFlattening.C\(11\)
    RUN={$MNRNF} VER=$MNDBVERSION SUFX=$MNDB root -q DoFlattening.C\(10\)
}

RUNNUMBER1=(${MNRNF})

function corr(){    ## only the first run
    RUN=${RUNNUMBER1[0]} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root DoRPAnalysis.C\($MXEVT\)
}


function allcorr(){ ## Go through from the second to the last
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	RUN=${RUN} VER=$MNVERSION SUFX=$MNDB DBVER=$MNDBVERSION root -b -q DoRPAnalysis.C >& $LOG 
	let I++
    done
}

echo ${RUNNUMBER1[0]}
grep function DoRPAnalysis.sh
env | grep MN

echo "type docorrection"

function docorrection() {
    flattening
    allcorr
}
