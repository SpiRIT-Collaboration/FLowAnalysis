RUN0={$RN132}#!/bin/bash                                                                                                                                           
source runList.sh

export DB=_rf.v5.1.0.cv0

#RUN0={$RNF112} DB0=$DB root calcFlw.C
#RUN0={$RN132s} DB0=$DB132 valgrind root calcFlw.C


RUN0={$RNF132} RUN1={$RNF108} RUN2={$RNF124} RUN3={$RNF112} DB0=$DB root calcFlw.C

