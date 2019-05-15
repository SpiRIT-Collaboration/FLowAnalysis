#! /bin/bash
#
# DoRPAnalysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh

source runList.sh
RUNNUMBER1=(${RNF132p})
RUNNUMBER1=(${RNF132r})
RUNNUMBER1=(${RNF132})
RUNNUMBER1=(${RNF108})
RUNNUMBER1=(${RNF124})
RUNNUMBER1=(${RNF112})

DBVERSION=25
VERSION=25.0
SUFX=BTt

##RUN=2841 SUFX=BTt VER=18.0 DBVER=18 root DoFlow_Analysis.C

function flowcorr(){
    RUN=${RUNNUMBER1[0]} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root DoRPAnalysis.C\($MXEVT\)
}


function allflowcorr(){
    typeset -i I=1
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	echo RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	echo RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	echo RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C : $LOG
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoRPAnalysis.C >& $LOG 
	let I++
    done
}