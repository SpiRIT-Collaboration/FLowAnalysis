#!/bin/bash

source setup.sh

#export MNMACRO=DoFlow_mod.C        ##<--- MACRO name
export MNMACRO=DoFlow_adv.C        ##<--- MACRO name

##---->>> EDIT here



##<----------- data
export MNRNF=$RNF108               ##<--- 
export MNDB=BTt                    ##<---
export MNVERSION=43.0              ##   <------@@ input 
export MNDBVERSION=$MNVERSION
##<----
export RPBS=0
export LC=
#20
export UC=
#40
export MNOVER=0
##------

##------------- RPSim  //Simulation
export MNDB=rpsim                  ##<---
##////----------
export MNVERSION=53                ##   <------@@ input 
export MNDBVERSION=$MNVERSION
export MNRNF=0403                  ##<--- 
##<-----------
export MNVERSION=50.7                ##   <------@@ input 
export MNDBVERSION=$MNVERSION
export MNRNF=0400                  ##<--- 
##////----------
export MNVERSION=22.0                ##   <------@@ input 
export MNDBVERSION=$MNVERSION
export MNRNF=0022                  ##<--- 
##---

echo $MNDBVERSION

export SUFX=$MNDB
export MNrunOne='SUFX=$MNDB  VER=$MNVERSION root $MACRO.C'

function dorpres()
{
#   RPBS=0 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNVERSION root DoRPResSubEvent.C\($1\)
   RPBS=0 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNVERSION root DoRPRes.C\($1\)
}


function doflowmulti()
{
    if [ -n "$1" ]; then
	export MNOVER=$1
	echo "output version -> "  $MNOVER
    elif [ -n "$2" ]; then
	export MNDBVERSION=$2
    fi

    PARTICLES=("2" "3" "4" "5" "6" "8")
    ##PARTICLES=("3")
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

    RPBS=0 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNOVER DBVER=$MNDBVERSION root -b -q $MNMACRO\($1\)
}

function doflow() 
{
    if [ -n "$2" ]; then
	export MNOVER=$2
	echo "output version -> "  $MNOVER
    elif [ -n "$3" ]; then
	export MNDBVERSION=$3
    fi

#    LC=0 UC=500 RPBS=$3 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNOVER root $MNMACRO\($1\)
#    LC=0 UC=80 RPBS=$3 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNOVER root $MNMACRO\($1\)
#    LC=0 UC=40 RPBS=$3 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNOVER root $MNMACRO\($1\)
    RPBS=0 RUN={$MNRNF} VER=$MNVERSION OUTVER=$MNOVER DBVER=$MNDBVERSION root $MNMACRO\($1\)
}

echo $MNVERSION to $MNOVER
cat DoFlow.sh |grep function
env|grep MN

echo "dorpres 0   :: Get centrality and Psi dependent correction"
echo "dorpres 1   :: Get Psi dependent correction parameter"
echo "dorpres 2   :: Get overall correction factor"
echo "doflow      :: open files "
echo "doflow 2  0(output version#) 1(Corretion versoin):: DoFlow_adv.C"
echo "-1:pi- 1:pi+, 2:p, 3:d, 4:t, 5:3He, 6: 4He, 7:n, 8:H" 
echo "Type  run #(partid) #(Output version)"


function doflowmdependence() 
{
    
    LC=0  UC=30  RUN=$MNRNF VER=$MNVERSION root -b -q $MNMACRO\(-4\)  
    LC=30 UC=40  RUN=$MNRNF VER=$MNVERSION root -b -q $MNMACRO\(-4\)
    LC=40 UC=50  RUN=$MNRNF VER=$MNVERSION root -b -q $MNMACRO\(-4\)
    LC=50 UC=80  RUN=$MNRNF VER=$MNVERSION root -b -q $MNMACRO\(-4\)

    LC=0  UC=30  RUN=$MNRNF VER=$MNVERSION OUTVER=5 root -b -q $MNMACRO\(8\)   
    LC=30 UC=40  RUN=$MNRNF VER=$MNVERSION OUTVER=6 root -b -q $MNMACRO\(8\) 
    LC=40 UC=50  RUN=$MNRNF VER=$MNVERSION OUTVER=7 root -b -q $MNMACRO\(8\) 
    LC=50 UC=80  RUN=$MNRNF VER=$MNVERSION OUTVER=8 root -b -q $MNMACRO\(8\) 
}
