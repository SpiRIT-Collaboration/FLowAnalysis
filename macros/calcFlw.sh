RUN0={$RN132}#!/bin/bash                                                                                                                                           
source runList.sh
source setup.sh

#export DB=_rf.v8.0.0.cv0
#export DB=_f0.v11
export DB=_rf.v11.0.1.cv2
export DB=_rf.v12.0.0.cv0
export DB=_f0.v14


#RUN0={$RN132} DB0=$DB132 valgrind root calcFlw.C
#RNF132="2841,2843,2844,2845,2846,2848,2849,2850,2851,2852,2855,2856,2857,2858,2859,2860,2861,2875,2877,2878,2879,2880,2881,2882,2883,2884,2887,2888,2889,2890,2891,2892,2893,2894,2896,2898,2899,2900,2901,2902,2903,2904,2905,2907,2914,2916,2917,2919,2920,2921,2922,2924,2925,2926,2927,2929,2930,2931,2932,2933,2934,2935,2936,2939,2940,2941,2942,2943,2944,2945,2946,2948,2955,2956,2958,2959,2960,2961,2962,2964,2965,2966,2968,2969,2970,2971,2972,2973,2975,2976,2977,2978,2979,2980,2981,2982,2983,2984,2985,2986,2988,2989,2990,2991"

echo 'type like this : '
export runFull='RUN0={$RNF132} RUN1={$RNF108} RUN2={$RNF124} RUN3={$RNF112} DB0=$DB root calcFlw.C'
export runOne='RUN0={"2841"} RUN1={"2273"} RUN2={"3061"} RUN3={"2542"} DB0=$DB root calcFlw.C'

echo ${runFull}
echo ${runOne} 
alias run132="RUN0={\$RNF132} DB0=$DB root calcFlw.C"
alias run108="RUN0={\$RNF108} DB0=$DB root calcFlw.C"
alias run124="RUN0={\$RNF124} DB0=$DB root calcFlw.C"
alias run112="RUN0={\$RNF112} DB0=$DB root calcFlw.C"
alias run132b="RUN0={\$RNF132} DB0=_rf.v10.2.5.cv5 root drawBootStrap.C"
alias run108s="RUN0={\$RNF108s} DB0=$DB root calcFlw.C"
alias run132s="RUN0={\$RNF132s} DB0=$DB root calcFlw.C"
alias run132n="UC=0 LC=1 RUN0={\$RNF132} DB0=$DB root calcFlw.C\(6\)"
alias |grep run1


export SPID=2
function exec132() {
    UC=0 LC=3 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=0 LC=1 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

function exec108() {
    UC=0 LC=3 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=0 LC=1 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

function exec124() {
    UC=0 LC=3 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=0 LC=1 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

function exec112() {
    UC=0 LC=3 RUN0={$RNF112} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=0 LC=1 RUN0={$RNF112} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF112} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF112} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

export SPID=6
function exec4() {
    UC=1 LC=2 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF112} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

export SPID=7
function azum() {
    AZ=1  RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
    AZ=10 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
    AZ=11 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
#   AZ=2  RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
#   AZ=20 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
#   AZ=21 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
#   AZ=3  RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
#   AZ=30 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
#   AZ=31 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
}


function azum12() {
    AZ=13  RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
    AZ=13  RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\)
}



cat calcFlw.sh |grep function


#PARTICLE=("2")
PARTICLE=("3" "4" "6" "7")
#PARTICLE=("5")
function multipleexec(){
    typeset -i I=0;
    while(( $I < ${#PARTICLE[@]} ))
    do
	SPID=${PARTICLE[I]}
	echo $I and $SPID 
	azum12
	let I++;
	if [ $I -ge ${#PARTICLE[@]} ]; then
	    break;
	fi
    done
}