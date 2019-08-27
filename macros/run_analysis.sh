#! /bin/bash
#
# run_analysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh

# TPC data
#
#VERSION=38
#export TPCDIR=/data/spdaq01/recoData/20190725/data/
#export RCVER=HEAD.1834.fc633ca

#VERSION=37
#export TPCDIR=/data/spdaq01/recoData/20190206/data/
#export RCVER=HEAD.1780.1e193e6

VERSION=39
export TPCDIR=/data/spdaq01/recoData/20190804/data/
export RCVER=HEAD.1853.e498ace

export ST132DIR=${STTPCDIR}
export ST108DIR=${STTPCDIR}
export ST124DIR=${STTPCDIR}
export ST112DIR=${STTPCDIR}
export STCkT100=${STTPCDIR}
export STCkT300=${STTPCDIR}

export STBBFITTER=db/BBFitter.root

#BigRIPS data
export STBEAM132=/cache/scr/spirit/DataAskedByMizuki/beam.Sn132_all/
export STBEAM108=/cache/scr/spirit/DataAskedByMizuki/beam.Sn108/
export STBEAM124=/cache/scr/spirit/DataAskedByMizuki/beam.Sn124/
export STBEAM112=/cache/scr/spirit/DataAskedByMizuki/beam.Sn112/

#KATANA data
export STKATANADIR=/xrootd/spdaq02/katana/root/katana/

#Kyoto data
export STKYOTODIR=/cache/scr/spirit/kaneko/rootfile/kyoto/
export STKYMLTDIR=/cache/scr/spirit/kaneko/rootfile/kyoto_re/mult/

#NeuLAND data
export STNLDIR=/cache/scr/spirit/NeuLand/neuland_4sep2018

#Anlaysis Flag
export STPC=1;
export BIGRIPS=1;
export KYOTOARRY=0;
export KATANA=1;
export NEULAND=1;

#--------------

source runList.sh

##--- for Process1 ------------------------------------
# *****> <Edit Here>
# Set RUNNUMBER1 
DBVERSION=0

RUNNUMBER1=(${RNF132})
#RUNNUMBER1=(${RNF108})
#RUNNUMBER1=(${RNF124})
#RUNNUMBER1=(${RNF112})
#RUNNUMBER1=(${RNFTEMP})
#RUNNUMBER1=(${RNF132r})
#RUNNUMBER1=(${RNF132s})

function execa() { ## Job for the 2841 with maximum event number =MEVT
    RUN=${RUNNUMBER1[0]} VER=$VERSION TPCDIR=$TPCDIR MXEVT=$MEVT DBVER=$DBVERSION root run_analysis.C 
}

MEVT=
function execb() {  ## batch job for the first run with maxium event number = MEVT
    echo $MEVT
    LOG=log/p1_${RUN}_v${VERSION}.log
    RUN=${RUNNUMBER1[0]} VER=$VERSION TPCDIR=$TPCDIR MXEVT=$MEVT DBVER=$DBVERSION root -b -q run_analysis.C >& $LOG &
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

#RUN=${RUNNUMBER1[0]}
#RUN=${RUN} VER=$VERSION root -b -q run_analysis.C\(1000\)

grep function run_analysis.sh

#execb
#execc







