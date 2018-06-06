RUN0={$RN132}#!/bin/bash                                                                                                                                           
source runList.sh

export DB=_rf.v5.1.0.cv0

RUN0={$RNF112} DB0=$DB root calcFlw.C
#RUN0={$RN132s} DB0=$DB132 valgrind root calcFlw.C


#RUN0={$RN132} RUN1={$RN108} RUN2={$RN124} RUN3={$RN112} DB0=$DB root calcFlw.C

