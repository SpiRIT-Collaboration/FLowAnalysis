RUN0={$RN132}#!/bin/bash                                                                                                                                           
source runList.sh
source setup.sh

#export DB=_rf.v8.0.0.cv0
#export DB=_f0.v11
export DB=_rf.v11.0
export DB=_rf.v11.0.1.cv2


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
alias run132b="RUN0={\$RNF132} DB0=_rf.v10.2.5.cv5 root drawBootStrap.C"
alias run108s="RUN0={\$RNF108s} DB0=$DB root calcFlw.C"
alias run132s="RUN0={\$RNF132s} DB0=$DB root calcFlw.C"
alias run132n="UC=0 LC=1 RUN0={\$RNF132} DB0=$DB root calcFlw.C\(6\)"
alias |grep run1


export SPID=2
function exec132() {
    UC=0 LC=1 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

function exec108() {
    UC=0 LC=1 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

function exec124() {
    UC=0 LC=1 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=1 LC=2 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
    UC=2 LC=3 RUN0={$RNF124} DB0=$DB root -b -q calcFlw.C\($SPID\) 
}

function exec112() {
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
cat calcFlw.sh |grep function


#PARTICLE=("2" "3" "4" "5")
PARTICLE=("2" "3" "5")
function multipleexec(){
    typeset -i I=0;
    while(( $I < ${#PARTICLE[@]} ))
    do
	SPID=${PARTICLE[I]}
	echo $I and $SPID 
#	UC=0 LC=3 RUN0={$RNF132} DB0=$DB root -b -q calcFlw.C\($SPID\)
	UC=0 LC=3 RUN0={$RNF108} DB0=$DB root -b -q calcFlw.C\($SPID\)
	let I++;
	if [ $I -ge ${#PARTICLE[@]} ]; then
	    break;
	fi
    done
}