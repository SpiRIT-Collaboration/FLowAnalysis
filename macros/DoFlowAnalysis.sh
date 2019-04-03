#! bin/bash

source ../build/config.sh

source runList.sh
RUNNUMBER1=(${RNF132})

DBVERSION=19
VERSION=19.1
SUFX=BTt

##RUN=2841 SUFX=BTt VER=18.0 DBVER=18 root DoFlow_Analysis.C

function flowcorr() {
    RUN=${RUNNUMBER1[0]} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root DoFlowAnalysis.C
}


function allflowcorr() {
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	RUN=${RUN} VER=$VERSION SUFX=$SUFX DBVER=$DBVERSION root -b -q DoFlowAnalysis.C
	let I++
	if ( $I -ge ${#RUNNUMBER1[@]} ); then
	    break;
	fi
    done
}
