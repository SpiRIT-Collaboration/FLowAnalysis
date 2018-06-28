RUN0={$RN132}#!/bin/bash                                                                                                                                           
source runList.sh

export DB=_rf.v5.1.0.cv0

#RUN0={$RN132s} DB0=$DB132 valgrind root calcFlw.C

echo 'type like this : '
export runFull='RUN0={$RNF132} RUN1={$RNF108} RUN2={$RNF124} RUN3={$RNF112} DB0=$DB root calcFlw.C'
export runOne='RUN0={"2841"} RUN1={"2273"} RUN2={"3061"} RUN3={"2542"} DB0=$DB root calcFlw.C'

echo ${runFull}
echo ${runOne} 
echo "RUN0={\$RNF132} DB0=$DB root calcFlw.C"
