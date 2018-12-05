#!/bin/bash                                                                                                                                           

source runList.sh

export DBP0=_rf.v11.0
export VER=2

#RUN0={$RNF08} DB0=$DBP2 VER=0 root  openFlw.C 
#RUN0={2841,$RNF132s,$RNF132r} DB0=$DBP0 VER=4 root  calibFlw.C 
#RUN0={$RNF132} DB0=$DBP0 VER=$VER root  calibFlw.C 
#RUN0={$RNF108} DB0=$DBP0 VER=$VER root  calibFlw.C 
RUN0={$RNF124} DB0=$DBP0 VER=$VER root -b -q calibFlw.C 
RUN0={$RNF112} DB0=$DBP0 VER=$VER root -b -q calibFlw.C 



