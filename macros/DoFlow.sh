source setup.sh

MNDB=BTt
MNVERSION=29.1
export MNSPID=2
export MNMACRO=DoFlow_adv.C

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
PARTICLE=("2" "3" "4" "5" "8")
export MNOVER=24
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

echo $MNVERSION to $MNOVER
cat DoFlow.sh |grep function
env|grep MN

