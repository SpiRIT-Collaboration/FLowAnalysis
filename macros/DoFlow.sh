#!/bin/bash

source setup.sh


RNF132="2841,2843,2844,2845,2846,2848,2849,2850,2851,2852,2855,2856,2857,2858,2859,2860,2861,2875,2877,2878,2879,2880,2881,2882,2883,2884,2887,2888,2889,2890,2891,2892,2893,2894,2896,2898,2899,2900,2901,2902"
#2903,2904,2905,2907,2914,2916,2917,2919,2920,2921,2922,2924,2925,2926,2927,2929,2930,2931,2932,2933,2934,2935,2936,2939,2940,2941,2942,2943,2944,2945,2946,2948,2955,2956,2958,2959,2960,2961,2962,2964,2965,2966,2968,2969,2970,2971,2972,2973,2975,2976,2977,2978,2979,2980,2981,2982,2983,2984,2985,2986,2988,2989,2990,2991,2992,2993,2997,2999,3000,3002,3003,3007,3039"

##---->>> EDIT here
export MNMACRO=DoFlow_adv.C        ##<--- MACRO name
export MNRNF={$RNF132}             ##<--- 
export MNDB=BTt                    ##<---
export MNVERSION=41.0              ##   <------@@ input 
export MNOVER=9
##<-----------

export MNrunOne='SUFX=$MNDB  VER=$MNVERSION root $MACRO.C'

function doflowmulti()
{
    if [ -n "$1" ]; then
	export MNOVER=$1
    fi
    echo Output Version is  $MNOVER

    PARTICLES=("3" "4" "5" "6" "8")
    typeset -i I=0;
    while(( $I < ${#PARTICLES[@]} ))
    do
	echo $I 
	doflowbatch ${PARTICLES[I]}
	let I++;
	if [ $I -ge ${#PARTICLES[@]} ]; then
	    break;
	fi
    done
}

function doflowbatch() 
{
    LC=0 UC=60 RUN=$MNRNF VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($1\)
}

function doflow() 
{
    if [ -n "$2" ]; then
	export MNOVER=$2
	echo $2 and $MNOVER
    fi

    LC=0 UC=60 RUN=$MNRNF VER=$MNVERSION OUTVER=$MNOVER root $MNMACRO\($1\)
}

echo $MNVERSION to $MNOVER
cat DoFlow.sh |grep function
env|grep MN
echo "doflow -2   :: Get centrality and Psi dependent correction"
echo "doflow      :: open files "
echo "doflow 2 0(output version#) :: DoFlow_adv.C"
echo "-1:pi- 1:pi+, 2:p, 3:d, 4:t, 5:3He, 6: 4He, 7:n, 8:H" 
echo "Type  run #(partid) #(Output version)"


function doflowmdependence() 
{
    LC=0  UC=45  RUN=$MNRNF VER=$MNVERSION OUTVER=4 root -b -q $MNMACRO\(8\) 
    LC=45 UC=55  RUN=$MNRNF VER=$MNVERSION OUTVER=5 root -b -q $MNMACRO\(8\) 
    LC=55 UC=65  RUN=$MNRNF VER=$MNVERSION OUTVER=6 root -b -q $MNMACRO\(8\) 
    LC=65 UC=100 RUN=$MNRNF VER=$MNVERSION OUTVER=7 root -b -q $MNMACRO\(8\)
}
