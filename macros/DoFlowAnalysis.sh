#! /bin/bash
#
# DoFlowAnalysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh

source runList.sh
RUNNUMBER1=(${RNF132})

DBVERSION=20
VERSION=20.1
SUFX=BTt

##RUN=2841 SUFX=BTt VER=18.0 DBVER=18 root DoFlow_Analysis.C

function flowcorr() {
    RUN=${RUNNUMBER1[0]} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root DoFlowAnalysis.C\($MXEVT\)
}


function allflowcorr() {
    typeset -i I=1
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoFlowAnalysis.C >& $LOG &
	let I++
	if ( $I -ge ${#RUNNUMBER1[@]} ); then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoFlowAnalysis.C >& $LOG &
	let I++
	if ( $I -ge ${#RUNNUMBER1[@]} ); then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoFlowAnalysis.C >& $LOG &
	let I++
	if ( $I -ge ${#RUNNUMBER1[@]} ); then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/p2_${RUN}_v${VERSION}.log
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoFlowAnalysis.C >& $LOG 
	let I++
    done
}
