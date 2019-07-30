source setup.sh

##---->>> EDIT here
export MNMACRO=DoFlow_adv.C        ##<--- MACRO name
export MNRNF={$RNF108}            ##<--- 
export MNDB=BTt                    ##<---
export MNVERSION=37.1              ##   <------@@ input 
export MNOVER=2
##<-----------

export MNrunOne='SUFX=$MNDB  VER=$MNVERSION root $MACRO.C'

function runmulti()
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

function runbatch() 
{
    LC=0 UC=55 RUN=$MNRNF VER=$MNVERSION OUTVER=$MNOVER root -b -q $MNMACRO\($1\)
}

function run() 
{
    LC=0 UC=55 RUN=$MNRNF VER=$MNVERSION OUTVER=$MNOVER root $MNMACRO\($1\)
}

echo $MNVERSION to $MNOVER
cat DoFlow.sh |grep function
env|grep MN
echo "At first, Do it run -1 "
