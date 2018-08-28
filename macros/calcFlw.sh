RUN0={$RN132}#!/bin/bash                                                                                                                                           
source runList.sh
source setup.sh

export DB=_rf.v8.0.0.cv0

#RUN0={$RN132} DB0=$DB132 valgrind root calcFlw.C

echo 'type like this : '
export runFull='RUN0={$RNF132} RUN1={$RNF108} RUN2={$RNF124} RUN3={$RNF112} DB0=$DB root calcFlw.C'
export runOne='RUN0={"2841"} RUN1={"2273"} RUN2={"3061"} RUN3={"2542"} DB0=$DB root calcFlw.C'

echo ${runFull}
echo ${runOne} 
alias run132="RUN0={\$RNF132} DB0=$DB root calcFlw.C"
alias run108="RUN0={\$RNF108} DB0=$DB root calcFlw.C"
alias run124="RUN0={\$RNF124} DB0=$DB root calcFlw.C"
alias run112="RUN0={\$RNF112} DB0=$DB root calcFlw.C"
alias run132s="RUN0={\$RNF132s} DB0=$DB root calcFlw.C"
alias |grep run1