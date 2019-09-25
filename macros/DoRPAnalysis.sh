#! /bin/bash
#
# DoRPAnalysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh
source runList.sh

#RNF132="2841"

##------>>> EDIT HERE 
export MNRNF=$RNF132
export MNINVERSION=41.2
export MNVERSION=41.2
export MNSFX=BTt
export REDO=0
##<----


function flattening() {
    RUN={$MNRNF} VER=$MNINVERSION DBVER=$MNVERSION SUFX=$MNSFX root -q -b DoFlattening.C\(11\)
    RUN={$MNRNF} VER=$MNINVERSION DBVER=$MNVERSION SUFX=$MNSFX root -q DoFlattening.C\(10\)
}

RUNNUMBER1=(${MNRNF})

function corr(){    ## only the first run
    RUN=${RUNNUMBER1[0]} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root DoRPAnalysis.C\($MXEVT\)
}


function allcorr(){ ## Go through from the second to the last
    typeset -i I=1
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} VER=$MNVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} VER=$MNVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	echo RUN=${RUN} VER=$MNVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${MNVERSION}.log
	RUN=${RUN} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C >& $LOG 
	let I++
    done
}

echo ${RUNNUMBER1[0]}
grep function DoRPAnalysis.sh
env | grep MN

echo "type flattenandcorrection <- all full process"
echo "type correction <- correction only"

function flattenandcorrection() {
    flattening
    RUN=${RUNNUMBER1[0]} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C
    allcorr
}


function correction() {
    RUN=${RUNNUMBER1[0]} VER=$MNINVERSION SUFX=$MNSFX DBVER=$MNVERSION root -b -q DoRPAnalysis.C
    allcorr
}
