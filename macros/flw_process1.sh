#! /bin/bash
#
# flw_process1.sh
# (c) Mizuki Nishimura
#
# ----------------------

# TPC data
#export STTPCDIR=/cache/scr/spirit/recoData/20180309/
#export STVERSION=1523.dc416ee   ###@develop.1535.1915a48
#export STTPCDIR=/cache/scr/spirit/recoData/20180719/
#export STVERSION=1665.2e3712e

export STTPCDIR=/xrootd/spdaq01/recoData/20180826/Sn132_JustinPRF_pion/
export STVERSION=GenieTemp.1728.18d377e
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
export STKATANADIR=/data/spdaq01/katana/root/katana/

#Kyoto data
export STKYOTODIR=/cache/scr/spirit/kaneko/rootfile/kyoto/
export STKYMLTDIR=/cache/scr/spirit/kaneko/rootfile/kyoto_re/mult/

#NeuLAND data
export STNLDIR=/cache/scr/spirit/NeuLand/neuland_23jul2018
#export STNLDIR=/cache/scr/spirit/NeuLand/neuland_3jul2018
#export STNLDIR=/cache/scr/spirit/NeuLand/neuland_18jun2018

#Anlaysis Flag
export STPC=1;
export BIGRIPS=1;
export KYOTOARRY=0;
export KATANA=0;
export NEULAND=1;

#--------------

#RUN=2900 VER=0 root -b -q flw_process1.C
#RUN=2841 VER=6 root  flw_process1.C

source runList.sh

##--- for Process1 ------------------------------------
# *****> <Edit Here>
# Set RUNNUMBER1 

#RUNNUMBER1=(${RNF132})
#RUNNUMBER1=(${RNF108})
#RUNNUMBER1=(${RNF124})
#RUNNUMBER1=(${RNF112})

#RUNNUMBER1=(${RNF132},${RNF108},${RNF124},${RNF112}) 
RUNNUMBER1="2900"

# *****> <Edit Here>     
VERSION=9

source setup.sh

function process1(){
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]} 
	LOG=log/prc1_${RUN}_v${VERSION}.log
	echo RUN=${RUN} VER=$VERSION root -b -q flw_process1.C '>&' $LOG
	RUN=${RUN} VER=$VERSION root -b -q flw_process1.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi

	RUN=${RUNNUMBER1[I]} 
	LOG=log/prc1_${RUN}_v${VERSION}.log
        echo RUN=${RUN} VER=$VERSION root -b -q flw_process1.C '>&' $LOG
	RUN=${RUN} VER=$VERSION root -b -q flw_process1.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi

	RUN=${RUNNUMBER1[I]} 
	LOG=log/prc1_${RUN}_v${VERSION}.log
        echo RUN=${RUN} VER=$VERSION root -b -q flw_process1.C '>&' $LOG
	RUN=${RUN} VER=$VERSION root -b -q flw_process1.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
            break;
        fi

	RUN=${RUNNUMBER1[I]} 
	LOG=log/prc1_${RUN}_v${VERSION}.log
        echo RUN=${RUN} VER=$VERSION root -b -q flw_process1.C '>&' $LOG
	RUN=${RUN} VER=$VERSION root -b -q flw_process1.C >& $LOG 
	let I++
    done
}


process1 





