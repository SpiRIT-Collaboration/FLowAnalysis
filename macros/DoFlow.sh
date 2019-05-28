source runList.sh
source setup.sh

DB=BTt
VERSION=25.0
export SPID=2
export runOne='RUN={"2841"} SUFX=$DB  VER=$VERSION root DoFlow.C'


function exec132() 
{
    UC=0 LC=3 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlow.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlow.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlow.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlow.C\($SPID\) 
}

function exec108() 
{
    UC=0 LC=3 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
}

function exec124() 
{
    UC=0 LC=3 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
}

function exec112() 
{
    UC=0 LC=3 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlow.C\($SPID\) 
}

export SPID=8
function exec4() 
{
#    UC=0 LC=80 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=8 root -b -q DoFlow.C\($SPID\) 
    LC=0 UC=80 RUN={$RNF108} SUFX=$DB  VER=$VERSION OUTVER=8 root -b -q DoFlow.C\($SPID\) 
    LC=0 UC=80 RUN={$RNF124} SUFX=$DB  VER=$VERSION OUTVER=8 root -b -q DoFlow.C\($SPID\) 
    LC=0 UC=80 RUN={$RNF112} SUFX=$DB  VER=$VERSION OUTVER=8 root -b -q DoFlow.C\($SPID\) 
}

export RNF={$RNF112}
export SPID=8
exprot OVER=25.0.2
function runonebyonebatch() 
{
    RUN=$RNF SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlow.C\($SPID\)
}

function runonebyone() 
{
    RUN=$RNF SUFX=$DB  VER=$VERSION OUTVER=$OVER root DoFlow.C\($SPID\)
}

function runonebyoneshort() 
{
    RUN={$RNF132r} SUFX=$DB  VER=$VERSION root DoFlow.C\($SPID\)
}



PARTICLE=("2" "3" "4" "6" "5") #"7")
function multipleexec()
{
    typeset -i I=0;
    while(( $I < ${#PARTICLE[@]} ))
    do
	SPID=${PARTICLE[I]}
	echo $I and $SPID 
	runbatch
	let I++;
	if [ $I -ge ${#PARTICLE[@]} ]; then
	    break;
	fi
    done
}

function runbatch() 
{
    RUN=$RNF SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlow.C\($SPID\)
}

cat DoFlow.sh |grep function

