RUN0={$RN132}#!/bin/bash                                                                                                                                           
source openRun.sh


export DB132=_rf.v4.0.5.cv3
export DB132=_rf.v4.5
export DB132=_rf.v4.6.0.cv0
export DB132=_mf.v4.7.0.cv0
export DB132=_rf.v4.7.0.cv0
export DB=_rf.v5.0.0.cv0
export DB=_rf.v5.1.0.cv0

RUN0={$RN132} DB0=$DB root calcFlw.C
#RUN0={$RN132s} DB0=$DB132 valgrind root calcFlw.C


#RUN0={$RN132} RUN1={$RN108} RUN2={$RN124} RUN3={$RN112} DB0=$DB root calcFlw.C

