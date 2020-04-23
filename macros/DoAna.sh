#! /BinRUN/Bash
#
# Dorpanalysis.Sh
# (C) Mizuki Nishimura
#
# ----------------------

source setup.sh
source runList.sh
export MNMACRO=DoFlow_adv.C          ##<--- Flow analysis MACRO name


##------>>> EDIT HERE 

##-- DATA
function data47()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNRNF=$RNF132
    export MNINVERSION=47            ##   <------@@ input
    export MNOUTVERSION=47.1         ##   <------@@ output 
    export MNOSUBV=
    export MNUC=40

    if [ -n "$1" ]; then
	SetStep3
    fi

    commonsetup
}
##<----
function data47.0()    ##CONFIG
{
    data47
    export MNOUTVERSION=47.0         ##   <------@@ output 
    commonsetup
}
##<----
function data49()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNRNF=$RNF132
    export MNINVERSION=49            ##   <------@@ input
    export MNOUTVERSION=49.0         ##   <------@@ output 
    export MNOSUBV=
    export MNUC=40

    if [ -n "$1" ]; then
	SetStep3
    fi

    commonsetup
##<----
} 
function data45()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNRNF=$RNF132
    export MNINVERSION=45            ##   <------@@ input
    export MNOUTVERSION=45.2         ##   <------@@ output 
    export MNOSUBV=
    export MNUC=40

    if [ -n "$1" ]; then
	SetStep3
    fi

    commonsetup
##<----
} 
function data44()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNRNF=$RNF132
    export MNINVERSION=44            ##   <------@@ input
    export MNOUTVERSION=44.0         ##   <------@@ output 
    export MNOSUBV=
    export MNUC=40

    if [ -n "$1" ]; then
	SetStep3
    fi

    commonsetup
##<----
}
function datapre()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNRNF=$RNF132
    export MNINVERSION=43            ##   <------@@ input
    export MNOUTVERSION=43.1         ##   <------@@ output 
    export MNOSUBV=
    export MNUC=40

    if [ -n "$1" ]; then
	SetStep3
    fi

    commonsetup
##<----
} ##----end DATA



##-- for simulation
function simtpc()   ##CONFIG
{
    export MNSFX=rpsim
    export MNUC=
    export MNRNF=0400 
    export MNINVERSION=50.7          ##   <------@@ input
    export MNOUTVERSION=50.8         ##   <------@@ output 
    if [ -n "$1" ]; then
	SetStep3
    fi
    commonsetup
}
function simfull()  ##CONFIG
{
    export MNSFX=rpsim
    export UC=
    export MNRNF=0022
    export MNINVERSION=22.0          ##   <------@@ input
    export MNOUTVERSION=22.5         ##   <------@@ output 
    if [ -n "$1" ]; then
	SetStep3
    fi
    commonsetup
}
##--end simulation


##-- common setup
function commonsetup()
{
    export IVER=$MNINVERSION
    export OVER=$MNOUTVERSION
    export SUFX=$MNSFX
    export DBVER=$MNINVERSION
    export OSUBV=$MNOSUBV
    export MXEVT=$MNMXEVT
    export REDO=$MNREDO
    export MNMXEVT=
    export MNREDO=
    export UC=$MNUC
    export RPBS=0

    RUNNUMBER1=(${MNRNF})
    anahelp
}
##<<< ---


##### analysis functions 
##>>>>>> Step 1 >
function flattening() {
    RUN={$MNRNF} root -q -b DoFlattening.C\(13\)
    RUN={$MNRNF} root -q -b DoFlattening.C\(12\)
    RUN={$MNRNF} root -q -b DoFlattening.C\(11\)
    RUN={$MNRNF} root -q DoFlattening.C\(10\)
}

##>>>>>> Step 2 >
function corr(){    ## The first run only; flatteing correction on interactive mode
    if [ $1 -eq 1 ]; then
	export REDO=1
	echo " Re-calculating Reaction plane #### "
    fi	

    RUN=${RUNNUMBER1[0]}  root DoRPAnalysis.C

    export REDO=
}

function corrb(){   ## The First run only; flattening correction on batch mode
    if [ $1 -eq 1 ]; then
	export REDO=1
	echo " Re-calculating Reaction plane #### "
    fi	

    RUN=${RUNNUMBER1[0]}  root -b -q DoRPAnalysis.C

    export REDO=
}


function allcorr(){ ## Go through from the second to the last
    if [ $1 -eq 1 ]; then
	export REDO=1
	echo " Re-calculating Reaction plane #### "
    fi	

    typeset -i I=1
    echo " Number of RUN : " ${#RUNNUMBER1[@]}
    while(( $I < ${#RUNNUMBER1[@]} ))
    do
	RUN=${RUNNUMBER1[I]}
	LOG=log/doana_${RUN}_v${MNOUTVERSION}.log
	echo $RUN 'root -b -q DoRPAnalysis.C  :' $LOG
	RUN=$RUN root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/doana_${RUN}_v${MNOUTVERSION}.log
	echo $RUN 'root -b -q DoRPAnalysis.C  :' $LOG
	RUN=$RUN root -b -q DoRPAnalysis.C >& $LOG &
	let I++
	if [ $I -ge ${#RUNNUMBER1[@]} ]; then
	    break;
	fi

	RUN=${RUNNUMBER1[I]}
	LOG=log/doana_${RUN}_v${MNOUTVERSION}.log
	echo $RUN 'root -b -q DoRPAnalysis.C  :' $LOG
	RUN=$RUN root -b -q DoRPAnalysis.C >& $LOG 
	let I++
    done

    export REDO=
}

function flattenandcorrection() {
    export REDO=0
    flattening
    corrb 
    allcorr
}

function correction() {
    corrb
    allcorr
}


##>>>>>> Step 3 >
function SetStep3()
{
    export RUN=${RUNNUMBER1[0]}  
    export IVER=$MNOUTVERSION
    export DBVER=$MNOUTVERSION
    showenv
}

function dorpres()
{
    SetStep3
    RUN={$MNRNF} root DoRPRes.C\($1\)
}

function doflowmdependence() 
{
    
    LC=0  UC=30  root -b -q $MNMACRO\(-4\)  
    LC=30 UC=40  root -b -q $MNMACRO\(-4\)
    LC=40 UC=50  root -b -q $MNMACRO\(-4\)
    LC=50 UC=80  root -b -q $MNMACRO\(-4\)

    LC=0  UC=30  OUTVER=5 root -b -q $MNMACRO\(8\)   
    LC=30 UC=40  OUTVER=6 root -b -q $MNMACRO\(8\) 
    LC=40 UC=50  OUTVER=7 root -b -q $MNMACRO\(8\) 
    LC=50 UC=80  OUTVER=8 root -b -q $MNMACRO\(8\) 
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
   RPBS=0 RUN={$MNRNF} DBVER=$MNOUTVERSION root -b -q $MNMACRO\($1\)
}

function doflow() 
{
    if [ -n "$2" ]; then
	export OSUBV=$2
    fi
    
    RPBS=0 RUN={$MNRNF} root $MNMACRO\($1\)
}

function doopen()
{
    RPBS=0 RUN={$MNRNF} root $MNMACRO\($1\)
}
##>>>>>> End of Step>



function anahelp()
{
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "**** Environment *********************** "
    echo "++++++++++++++++++++++++++++++++++++++++"
    env | grep MN
    showenv
    echo "RUN --------------> " ${RUNNUMBER1[0]}
    echo "Max event number --> \$MXEVT " $MXEVT
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "type anahelp()"
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo " Step 0 > corr 1 OR correction 1 :: Re-caluculate reaction plane "
    echo " Step 1 > flattening OR flattenandcorrection "
    echo " Step 2 > corr(the first run) OR allcorr(excpt the first run) OR  correction(full)"
    echo " Step 3 > dorpres 0 "
    echo "          dorpres 0   :: Get centrality and Psi dependent correction"
    echo "          dorpres 1   :: Get Psi dependent correction parameter"
    echo "          dorpres 2   :: Get overall correction factor"
    echo " Step 4 > doflow       :: open files "
    echo "          doflow 2  0(output version#) 1(Corretion versoin):: DoFlow_adv.C"
    echo "                  -1:pi- 1:pi+, 2:p, 3:d, 4:t, 5:3He, 6: 4He, 7:n, 8:H" 
    echo "                  run #(partid) #(Output version)"
}

function showenv()
{
    echo "----------------------------------------"
    echo run number    :: $MNRNF
    echo Sufix         :: $SUFX
    echo Input         :: $IVER
    echo Output        :: $OVER
    echo DataBase      :: $DBVER
    echo Output subVer :: $RPBS
    echo Max maltiplicity :: $UC
    echo Min multiplicity :: $LC
    echo "++++++++++++++++++++++++++++++++++++++++"
}

grep "CONFIG"  DoAna.sh 
anahelp

