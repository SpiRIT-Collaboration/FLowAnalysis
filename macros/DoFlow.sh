source setup.sh

MNDB=BTt                           ##<---
MNVERSION=35.0                     ##   <------@@ input 
export MNSPID=2                    ##<--- pid: 2:p, 3:d, 4:t, 5:3He, 6:4He, 7:n, 8:n
export MNMACRO=DoFlow_adv.C        ##<---

export MNrunOne='SUFX=$MNDB  VER=$MNVERSION root $MACRO.C'


function exec4systems()   ## 4 systems
{
    export MNSPID=8
    LC=0  UC=100 RUN={$RNF132} SUFX=$DB  VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($MNSPID\) 
    LC=0  UC=100 RUN={$RNF108} SUFX=$DB  VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($MNSPID\) 
    LC=0  UC=100 RUN={$RNF124} SUFX=$DB  VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($MNSPID\) 
    LC=0  UC=100 RUN={$RNF112} SUFX=$DB  VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($MNSPID\) 
}


function runonebyonebatch() 
{
    RUN=$MNRNF SUFX=$DB  VER=$MNVERSION OUTVER=$MNOVER root -b -q $MACRO\($MNSPID\)
}

function runonebyone() 
{
    RUN=$MNRNF SUFX=$DB  VER=$MNVERSION OUTVER=$MNOVER root $MACRO\($MNSPID\)
}

function runonebyoneshort() 
{
    RUN={$MNRNF132r} SUFX=$DB  VER=$MNVERSION root $MACRO\($MNSPID\)
}


export MNRNF={$RNF108}
export MNRNF={$RNF132}
#PARTICLE=("3" "4" "5" "8")
PARTICLE=("2")
#export MNOVER=0
function multipleexec()
{
    typeset -i I=0;
    while(( $I < ${#PARTICLE[@]} ))
    do
	MNSPID=${PARTICLE[I]}
	echo $I and $MNSPID 
	runbatch
	let I++;
	if [ $I -ge ${#PARTICLE[@]} ]; then
	    break;
	fi
    done
}

function runbatch() 
{
    LC=0 UC=55 RUN=$MNRNF VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($MNSPID\)
}

function run() 
{
    LC=0 UC=55 RUN=$MNRNF VER=$MNVERSION OUTVER=$MNOVER root $MNMACRO
}

echo $MNVERSION to $MNOVER
cat DoFlow.sh |grep function
env|grep MN

