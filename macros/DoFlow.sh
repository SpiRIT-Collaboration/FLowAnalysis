#!/bin/bash

source setup.sh

##---->>> EDIT here
export MNMACRO=DoFlow_adv.C        ##<--- MACRO name
export MNRNF={$RNF132}             ##<--- 
export MNDB=BTt                    ##<---
export MNVERSION=40.0              ##   <------@@ input 
export MNOVER=1
##<-----------

export MNrunOne='SUFX=$MNDB  VER=$MNVERSION root $MACRO.C'

function doflowmulti()
{
    PARTICLES=("3" "4" "5")
    typeset -i I=0;
    while(( $I < ${#PARTICLES[@]} ))
    do
	echo $I 
	runbatch ${PARTICLES[I]}
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
echo "doflow -1 :: Get centrality and Psi dependent correction"
echo "doflow    :: open files "
echo "doflow  2 0 :: DoFlow_adv.C"
echo "Type  run #(partid) #(Output version)"