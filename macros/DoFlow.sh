source runList.sh
source setup.sh

DB=BTt
VERSION=23.1
OVER=0
export SPID=2
export runOne='RUN={"2841"} SUFX=$DB  VER=$VERSION root DoFlowAnalysis.C'


function exec132() 
{
    UC=0 LC=3 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlowAnalysis.C\($SPID\) 
}

function exec108() 
{
    UC=0 LC=3 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
}

function exec124() 
{
    UC=0 LC=3 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
}

function exec112() 
{
    UC=0 LC=3 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=0 LC=1 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=2 LC=3 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
}

export SPID=2
function exec4() 
{
    UC=1 LC=2 RUN={$RNF132} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF108} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF124} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
    UC=1 LC=2 RUN={$RNF112} SUFX=$DB  VER=$VERSION root -b -q DoFlowAnalysis.C\($SPID\) 
}

function runonebyonebatch() 
{
    RUN={$RNF132} SUFX=$DB  VER=$VERSION OUTVER=$OVER root -b -q DoFlowAnalysis.C\($SPID\)
}

function runonebyone() 
{
    RUN={$RNF132} SUFX=$DB  VER=$VERSION root DoFlowAnalysis.C\($SPID\)
}

function runonebyoneshort() 
{
    RUN={$RNF132r} SUFX=$DB  VER=$VERSION root DoFlowAnalysis.C\($SPID\)
}


PARTICLE=("3" "4" "6" "5") #"7")
function multipleexec()
{
    typeset -i I=0;
    while(( $I < ${#PARTICLE[@]} ))
    do
	SPID=${PARTICLE[I]}
	echo $I and $SPID 
	runonebyonebatch
	let I++;
	if [ $I -ge ${#PARTICLE[@]} ]; then
	    break;
	fi
    done
}

cat DoFlowAnalysis.sh |grep function

