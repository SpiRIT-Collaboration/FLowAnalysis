#! /bin/bash
#
# flw_process1.sh
# (c) Mizuki Nishimura
#
# ----------------------

# TPC data
export STTPCDIR=/cache/scr/spirit/recoData/20180309/
export ST132DIR=${STTPCDIR}132Sn124Sn/
export ST108DIR=${STTPCDIR}108Sn112Sn/
export ST124DIR=${STTPCDIR}124Sn112Sn/
export ST112DIR=${STTPCDIR}112Sn124Sn/
export STCkT100=${STTPCDIR}Cocktail100MeV/
export STCkT300=${STTPCDIR}Cocktail300MeV/

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
export STNLDIR=/cache/scr/spirit/NeuLand/

#Anlaysis Flag
export BIGRIPS=1;
export KYOTOARRY=0;
export KATANA=0;
export NEULAND=1;

#--------------

#RUN=2900 VER=0 root -b -q flw_process1.C
#RUN=2841 VER=6 root  flw_process1.C

### look at procRun.sh

##--- for Process1 ------------------------------------
# *****> <Edit Here>
# Set run number you want to analize. 
# See procRun.sh
RUNNUMBER1=(\
"2542" "2543" "2544" "2545" "2546" "2547" "2548" "2549" \
"2550" "2551" "2552" "2553" "2554" "2555" "2556" "2557" "2558" "2559" \
"2560" "2562" "2563" "2564" "2565" "2566" "2567" "2568" "2569" "2570" \
"2571" "2572" "2573" "2574" "2575" "2576" "2578" "2579" "2580" "2581" \
"2582" "2583" "2584" "2585" "2586" "2587" "2588" "2589" "2590" "2591" \
"2592" "2593" "2594" "2595" "2596" "2597" "2598" "2599" "2600" "2601" \
"2603" "2604" "2605" "2607" "2609" "2610" \
"2613" "2617" "2618" "2619" "2620" "2621" "2622" "2623" \
"2624" "2632" "2633" "2634" "2643" "2645" "2646" "2647" "2648" \
"2649" "2650" "2652" "2653" \
)

VERSION=5


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





