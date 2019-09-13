#!/bin/bash

source setup.sh

##---->>> EDIT here
export MNMACRO=DoFlow_adv.C        ##<--- MACRO name
export MNRNF={$RNF108}             ##<--- 
export MNDB=BTt                    ##<---
export MNVERSION=40.0              ##   <------@@ input 
export MNOVER=1
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
