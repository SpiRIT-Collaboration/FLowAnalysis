#! /BinRUN/Bash
#
# Dorpanalysis.Sh
# (C) Mizuki Nishimura
#
# ----------------------

source setup.sh
source runList.sh
export MNMACRO=DoFlow_adv.C          ##<--- Flow analysis MACRO name

unset IVER
unset OVER
unset SUFX
unset DBVER
unset OSUBV
unset MXEVT
unset REDO
unset MNMXEVT
unset UC
unset RPBS
unset RUNNUMBER1
unset MNRNF
unset MNINVERSION
unset MNOUTVERSION

##------>>> EDIT HERE 
export MNUC=40
export MNLC=0

##-- DATA
function data132()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##-- 
    export MNSN=132
    export MNRNF=$RNF132
    export MNINVERSION=$1            ##   <------@@ input
    export MNOUTVERSION=$2           ##   <------@@ output 
    export MNOSUBV=

    commonsetup
##<----
}

function data108()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNSN=108
    export MNRNF=$RNF108
    export MNINVERSION=$1            ##   <------@@ input
    export MNOUTVERSION=$2           ##   <------@@ output 
    export MNOSUBV=

    commonsetup
##<----
}

function data112()    ##CONFIG
{
    export MNSFX=BTt
    export MDCUT= #0.0               ##   <------@@ mid-rapidity cut
##--
    export MNSN=112
    export MNRNF=$RNF112
    export MNINVERSION=$1            ##   <------@@ input
    export MNOUTVERSION=$2           ##   <------@@ output 
    export MNOSUBV=

    commonsetup
##<----
}

##-- for simulation
function simtpc()   ##CONFIG
{
    export MNSFX=rpsim
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
    export MNMXEVT=
    export UC=$MNUC
    export LC=$MNLC
    export RPBS=0

    RUNNUMBER1=(${MNRNF})
    showenv
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
    if [ -n "$1" ]; then
	export REDO=1
	echo " Re-calculating Reaction plane #### "
    fi	
   
    RUN=${RUNNUMBER1[0]}  root DoRPAnalysis.C

    export REDO=
}

function corrb(){   ## The First run only; flattening correction on batch mode
    if [ -n "$1" ]; then
	export REDO=1
	echo " Re-calculating Reaction plane #### "
    fi	

    RUN=${RUNNUMBER1[0]}  root -b -q DoRPAnalysis.C

    export REDO=
}


function allcorr(){ ## Go through from the second to the last
    if [ -n "$1" ]; then
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
    RUN={$MNRNF} root DoRPRes.C\($1\)
}

function doflowmdependence() 
{
#    LC=0   UC=30   RUN={$MNRNF} root -b -q DoRPRes.C\(0\)
#    LC=30  UC=40   RUN={$MNRNF} root -b -q DoRPRes.C\(0\)
#    LC=40  UC=50   RUN={$MNRNF} root -b -q DoRPRes.C\(0\)
#    LC=50  UC=80   RUN={$MNRNF} root -b -q DoRPRes.C\(0\)

    LC=0  UC=30 OSUBV=15 doflowbatch $1
    LC=30 UC=40 OSUBV=16 doflowbatch $1
    LC=40 UC=50 OSUBV=17 doflowbatch $1
    LC=50 UC=80 OSUBV=18 doflowbatch $1
}

function doflowmulti()
{
    if [ -n "$1" ]; then
	export OSUBV=$1
	echo "output version -> "  $MNOVER
    fi

    PARTICLES=("3" "4" "5" "6")
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
   RPBS=0 RUN={$MNRNF} root -b -q $MNMACRO\($1\)
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
    echo "----------------------------------------"
    echo "RUN --------------> " ${RUNNUMBER1[0]}
    echo "Max event number --> \$MXEVT " $MXEVT
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "type anahelp()"
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo " Setup >> data132/data108  ##.#(INverseion) ##.#(OutVersion)"
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
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "**** Environment *********************** "
    echo "++++++++++++++++++++++++++++++++++++++++"
    env | grep MN
    echo "----------------------------------------"
    echo BEAM          :: $MNSN
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

