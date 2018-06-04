#!/bin/bash                                                                                                                                           

export RNF132t=2843
export RNF132s=2841,2843,2844,2845,2846,2848,2849,2850,2851,2852
export RNF108s=2273,2274,2275,2276,2283,2284,2285,2286,2288,2289
export RNF124s=3059,3061,3062,3065,3066,3068,3069,3071,3074,3075

export RNm132=2841,2843,2844,2845,2846,2848,2849,2850,2851 #,2852,2855,2856,2857,2858,2859,2860,2861,2875


source openRun.sh

export DBP2=_rf.v5.1

RUN0={$RN108} DB0=$DBP2 VER=0 root  calibFlw.C 


#RUN0={$RN132t} DB0=$DBP2 VER=3 root -b -q calibFlw.C > logz


