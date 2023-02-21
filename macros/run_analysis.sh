#! /bin/bash
#
# run_analysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh

# TPC data
#
VERSION=3
export TPCDIR=/home/recoData/20200529/data
export RCVER=develop.1988.bf2b00e
export SUFX=flow
###export SUFX=ph

##g4 data
#export TPCDIR=/home/isobe/20200909MizukiGeant4/recoData.20200910/
#mizuki_000003_s0.reco.v1.04.root

##export STBBFITTER=db/BBFitter.root
export STBBFITTER=db/FlattenPID.root # Kaneko-kun's PID calibration files done with 20190804/data.corr_ExBl112/ 

#BigRIPS data
export STBEAM132=/home/spana01/DataAskedByMizuki/beam.Sn132_all/
export STBEAM108=/home/spana01/DataAskedByMizuki/beam.Sn108/
export STBEAM124=/home/spana01/DataAskedByMizuki/beam.Sn124/
export STBEAM112=/home/spana01/DataAskedByMizuki/beam.Sn112/

#KATANA data

export STKATANADIR=/xrootd/spdaq02/katana/root/katana/

#Kyoto data
export STKYOTODIR=/cache/scr/spirit/kaneko/rootfile/kyoto/
export STKYMLTDIR=/cache/scr/spirit/kaneko/rootfile/kyoto_re/mult/

#NeuLAND data
export STNLDIR=/cache/scr/spirit/NeuLand/neuland_4sep2018

#Anlaysis Flag
export STPC=1;
export BIGRIPS=0;
export KYOTOARRY=0;
export KATANA=0;
export NEULAND=0;

#--------------

source runList.sh
##--- for Process1 ------------------------------------
# *****> <Edit Here>
# Set RUNNUMBER1 
DBVERSION=0

#export TPCDIR=$TPCDIR

#RNF132="2894"
RUNNUMBER1=(${RNF132})
#RUNNUMBER1=(${RNF108})
#RUNNUMBER1=(${RNF124})
#RUNNUMBER1=(${RNF112})
#RUNNUMBER1=(${RNFTEMP})
#RUNNUMBER1=(${RNF132ss})
#RUNNUMBER1=(${RNF108ss})
#RUNNUMBER1=(${RNFG4})

#RUNNUMBER1="2595"

echo $TPCDIR

function makeflow(){   ## batch job from the second to the end.
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do

	RUN=${RUNNUMBER1[I]} 
	echo RUN=${RUN} VER=$VERSION root -b -q makeFlowBranch.C
	RUN=${RUN} VER=$VERSION root -b -q makeFlowBranch.C >& log/error_${RUN}.log
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi
    done
}
function makephys(){   ## batch job from the second to the end.
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do

	RUN=${RUNNUMBER1[I]} 
	echo RUN=${RUN} root -b -q makePhysicsTree.C
	RUN=${RUN} root -b -q makePhysicsTree.C >& log/error_${RUN}.log
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi
    done
}

function exe3(){   ## batch job from the second to the end.
    typeset -i I=2
    while(( $I < ${#RUNNUMBER1[@]} ))
    do

	RUN=${RUNNUMBER1[I]} 
	echo RUN=${RUN} root -b -q exe3.C
	RUN=${RUN} root -b -q exe3.C
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi
    done
}

function execa() { ## Job for the 2841 with maximum event number = MXEVT
    RUN=${RUNNUMBER1[0]} VER=$VERSION TPCDIR=$TPCDIR DBVER=$DBVERSION MXEVT=$1 root run_analysis.C 
}


function execb() {  ## batch job for the first run with maxium event number = MXEVT
    echo $MXEVT
    RUN=${RUNNUMBER1[0]} 
    LOG=log/p1_${RUN}_v${VERSION}.log
    RUN=${RUNNUMBER1[0]} VER=$VERSION TPCDIR=$TPCDIR  DBVER=$DBVERSION MXEVT=$1 root -b -q run_analysis.C >& $LOG &
}

function execc(){   ## batch job from the second to the end.
    typeset -i I=1
    while(( $I < ${#RUNNUMBER1[@]} ))
    do

	RUN=${RUNNUMBER1[I]} 
	LOG=log/p1_${RUN}_v${VERSION}.log
	echo RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C '>&' $LOG
	RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi

	RUN=${RUNNUMBER1[I]} 
	LOG=log/p1_${RUN}_v${VERSION}.log
        echo RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C '>&' $LOG
	RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi

	RUN=${RUNNUMBER1[I]} 
	LOG=log/p1_${RUN}_v${VERSION}.log
        echo RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C '>&' $LOG
	RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi

	RUN=${RUNNUMBER1[I]} 
	LOG=log/p1_${RUN}_v${VERSION}.log
        echo RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C '>&' $LOG
	RUN=${RUN} VER=$VERSION DBVER=$DBVERSION root -b -q run_analysis.C >& $LOG 
	let I++
    done
}

function execd() { ## Job for the 2841 with maximum event number = MXEVT
    RUN=2843 VER=$VERSION TPCDIR=$TPCDIR DBVER=$DBVERSION MXEVT=$1 root run_analysis.C 
    RUN=2898 VER=$VERSION TPCDIR=$TPCDIR DBVER=$DBVERSION MXEVT=$1 root run_analysis.C 
}

#RUN=${RUNNUMBER1[0]}
#RUN=${RUN} VER=$VERSION root -b -q run_analysis.C\(1000\)

grep function run_analysis.sh


###--- show envirounment
echo run$RUNNUMBER1 with ${#RUNNUMBER1[@]} runs 
echo v$VERSION
##----






