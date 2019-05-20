#! /bin/bash
#
# run_analysis.sh
# (c) Mizuki Nishimura
#
# ----------------------

source ../build/config.sh

# TPC data

export TPCDIR=/xrootd/spdaq02/recoData/20190206/data/
export RCVER=HEAD.1780.1e193e6
#export TPCDIR=/xrootd/spdaq02/recoData/20181219/data/
#export RCVER=HEAD.1769.ef17b59
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
DBVERSION=25
VERSION=27
MEVT=100



function execa() {
    RUNNUMBER1="2841"
    MXEVT=$MEVT
    RUN=${RUNNUMBER1[0]} VER=$VERSION TPCDIR=$TPCDIR MXEVT=$MXEVT DBVER=$DBVERSION root run_analysis.C 
}


RUNNUMBER1=(${RNF132})
#RUNNUMBER1=(${RNF108})
#RUNNUMBER1=(${RNF124})
#RUNNUMBER1=(${RNF112})
#RUNNUMBER1=(${RBF132} ${RNF108} ${RNF124} ${RNF112}) 
#RUNNUMBER1=(${RNFTEMP})
#RUNNUMBER1=(${RNF132r})
#RUNNUMBER1=(${RNF132p})

MXEVT=
function execb() {
    LOG=log/p1_${RUN}_v${VERSION}.log
    RUN=${RUNNUMBER1[0]} VER=$VERSION TPCDIR=$TPCDIR MXEVT=$MXEVT DBVER=$DBVERSION root -b -q run_analysis.C >& $LOG &
}

function exec(){
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







