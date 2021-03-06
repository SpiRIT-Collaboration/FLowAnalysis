#! /Bin/Bash
#
# Dorpanalysis.Sh
# (C) Mizuki Nishimura
#
# ----------------------

source setup.sh
source runList.sh
export MNMACRO=DoFlow_fin.C          ##<--- Flow analysis MACRO name

unset IVER
unset OVER
unset SUFX
unset DBVER
unset OSUBV
unset MXEVT
unset REDO
unset MNMXEVT
unset TRK
unset UC
unset LC
unset RPBS
unset RUNNUMBER1
unset MNRNF
unset MNINVERSION
unset MNOUTVERSION
unset CCPSI
unset PHICUT

##------>>> EDIT HERE 
#export MNLC=0
#export MNUC=50

#export MNLC=55
#export MNUC=80

##-->> b3fm
##-->(132:52> 2.9fm)
##-->(132:42> 5.2fm)

multL=" 0,20,35,40,50,55,42,42,40,42, 0, 0"
multU="35,40,45,55,65,80,52,56,54,55,54,55"
##      0  1  2  3  4  5  6  7  8  9 10 11

export MNTRK=3
export PHICUT=2 #phi>135
export PHICUT=1 #phi<45

#########----------------------
MLC=(${multL})
MUC=(${multU})
export MNLC=${MLC[$MNTRK]}
export MNUC=${MUC[$MNTRK]}


export MNSFX=BTt
export MDCUT=0. #0.0               ##   <------@@ mid-rapidity cut
export MNINVERSION=52.15           ##   <------@@ input
export MNOUTVERSION=52.15          ##   <------@@ output 


##-- DATA
function data108s()    ##CONFIG
{
##--
    export MNSN=108
    export MNRNF=$RNF108s
    export MNOSUBV=

    if [ -n "$1" ]; then
	export MNINVERSION=$1            ##   <------@@ input
	export MNOUTVERSION=$2           ##   <------@@ output 
    fi

    commonsetup
##<----
}

function data132()    ##CONFIG
{
##-- 
    export MNSN=132
    export MNRNF=$RNF132
    export MNMXEVT=
    export MNOSUBV=

    if [ -n "$1" ]; then
	export MNINVERSION=$1            ##   <------@@ input
	export MNOUTVERSION=$2           ##   <------@@ output 
    fi

    commonsetup
##<----
}

function data108()    ##CONFIG
{
##--
    export MNSN=108
    export MNRNF=$RNF108
    export MNOSUBV=

    if [ -n "$1" ]; then
	export MNINVERSION=$1            ##   <------@@ input
	export MNOUTVERSION=$2           ##   <------@@ output 
    fi

    commonsetup
##<----
}

function data124()    ##CONFIG
{
##--
    export MNSN=124
    export MNRNF=$RNF124
    export MNOSUBV=

    if [ -n "$1" ]; then
	export MNINVERSION=$1            ##   <------@@ input
	export MNOUTVERSION=$2           ##   <------@@ output 
    fi

    commonsetup
##<----
}

function data112()    ##CONFIG
{
##--
    export MNSN=112
    export MNRNF=$RNF112
    export MNOSUBV=

    if [ -n "$1" ]; then
	export MNINVERSION=$1            ##   <------@@ input
	export MNOUTVERSION=$2           ##   <------@@ output 
    fi

    commonsetup
##<----
}


##-- for simulation
function simtpc()   ##CONFIG
{
    export MNSN=100
    export MNSFX=rpsim
    export MNRNF=0022 
    export MNINVERSION=22.0          ##   <------@@ input
    export MNOUTVERSION=22.0         ##   <------@@ output 

    commonsetup
}
function simfull()  ##CONFIG
{
    export MNSN=100
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
#    export DBVER=$MNOUTVERSION
    export OSUBV=$MNOSUBV
    export MXEVT=$MNMXEVT
    export MNMXEVT=
    export UC=$MNUC
    export LC=$MNLC
    export RPBS=0
    export CCPSI=1 #phi dependent correction

    export RUN={$MNRNF} 
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
function allcorr() {
    corrb $1
    restcorr $1
}

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

function restcorr(){ ## Go through from the second to the last
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
    RUN={$MNRNF} root -q DoRPRes.C\($1\)

}

function doflowmid132() 
{
    data132
    PARTICLES=("6" "5" "4" "3" "2" )
    LC=45  UC=57 PHICUT=1 doflowmulti 122
    LC=45  UC=57 PHICUT=2 doflowmulti 123
}
function doflowmid108() 
{
    data108
    PARTICLES=("6" "5" "4" "3" "2" )
    LC=44  UC=55 PHICUT=1 doflowmulti 122
    LC=44  UC=55 PHICUT=2 doflowmulti 123
}
function doflowmid112() 
{
    data112
    PARTICLES=("6" "5" "4" "3" "2" )
    LC=44  UC=55 PHICUT=1 doflowmulti 122
    LC=44  UC=54 PHICUT=2 doflowmulti 123
}

function doflowmdependence() 
{
    PARTICLES=("6" "5" "4" "3" "2" )
    LC=0  UC=20 PHICUT=1 doflowmulti 100
    LC=0  UC=20 PHICUT=2 doflowmulti 101
    LC=20 UC=30 PHICUT=1 doflowmulti 102
    LC=20 UC=30 PHICUT=2 doflowmulti 103
    LC=30 UC=35 PHICUT=1 doflowmulti 104
    LC=30 UC=35 PHICUT=2 doflowmulti 105
    LC=35 UC=40 PHICUT=1 doflowmulti 106
    LC=35 UC=40 PHICUT=2 doflowmulti 107
    LC=40 UC=45 PHICUT=1 doflowmulti 108
    LC=40 UC=45 PHICUT=2 doflowmulti 109
    LC=44 UC=48 PHICUT=1 doflowmulti 110
    LC=44 UC=48 PHICUT=2 doflowmulti 111
    LC=45 UC=50 PHICUT=1 doflowmulti 112
    LC=45 UC=50 PHICUT=2 doflowmulti 113
    LC=50 UC=55 PHICUT=1 doflowmulti 114
    LC=50 UC=55 PHICUT=2 doflowmulti 115
    LC=55 UC=60 PHICUT=1 doflowmulti 116
    LC=55 UC=60 PHICUT=2 doflowmulti 117
    LC=60 UC=65 PHICUT=1 doflowmulti 118
    LC=60 UC=65 PHICUT=2 doflowmulti 119
    LC=65 UC=80 PHICUT=1 doflowmulti 120
    LC=65 UC=80 PHICUT=2 doflowmulti 121
}

function doflowmultiHeavy()
{
    PARTICLES=("3" "4" "5" "6")
    doflowmulti $1
}

function doflowmultiall()
{
    PARTICLES=("2" "3" "4" "5" "6")
    doflowmulti $1
}

function doflowmulti()
{
    if [ -n "$1" ]; then
	export OSUBV=$1
	echo "output version -> "  $OSUBV
    fi


#    PARTICLES=("2" "3" "4" "5" "6")
    typeset -i I=0;
    while(( $I < ${#PARTICLES[@]} ))
    do
	echo $I 
	doflowbatch ${PARTICLES[I]} $1 
	let I++;
	if [ $I -ge ${#PARTICLES[@]} ]; then
	    break;
	fi
    done
}


function doflowbatch() 
{
    if [ -n "$2" ]; then
	export OSUBV=$2
    fi
   RPBS=0 RUN={$MNRNF} root -b -q $MNMACRO\($1\) 
   unset OSUBV
}

function doflow() 
{
    if [ -n "$2" ]; then
	export OSUBV=$2
    fi
    
    RPBS=0 RUN={$MNRNF} root $MNMACRO\($1\)
    unset OSUBV
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
    echo " Step 0 > corr 1 OR allcorr 1 :: Re-caluculate reaction plane "
    echo " Step 1 > flattening OR flattenandcorrection "
    echo " Step 2 > corr(the first run) OR restcorr(excpt the first run) OR  allcorr(full)"
#    echo " Step 3 > dorpres 0 "
#    echo "          dorpres 0   :: Get centrality and Psi dependent correction"
#    echo "          dorpres 1   :: Get Psi dependent correction parameter"
#    echo "          dorpres 2   :: Get overall correction factor"
    echo " Step 3 > doflow       :: open files "
    echo "          doflow 2  0(output version#) 1(Corretion versoin)::" $MNMACRO
    echo "                  -1:pi- 1:pi+, 2:p, 3:d, 4:t, 5:3He, 6: 4He, 7:n, 8:H" 
    echo "                  run #(partid) #(Output version)"
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "Change runnumber : export MNRNF=\$RNF132s then commonsetup"
}

function showenv()
{
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "**** Environment *********************** "
    echo "++++++++++++++++++++++++++++++++++++++++"
    env | grep MN
    echo "----------------------------------------"
    echo "BEAM             ::"  $MNSN
    echo "run number       ::"  $MNRNF
    echo "Sufix            ::"  $SUFX
    echo "Input            ::"  $IVER
    echo "Output           ::"  $OVER
    echo "DataBase         ::"  $DBVER
    echo "Output subVer    ::"  $RPBS
    echo "Selected track   ::"  $TRK
    echo "Number of event  ::"  $MXEVT
    echo "Mid-Y cut        ::"  $MDCUT
    echo "Min multiplicity ::"  $LC
    echo "Max maltiplicity ::"  $UC
    echo "Phi Angle cut    ::"  $PHICUT
    echo "++++++++++++++++++++++++++++++++++++++++"
}

function showdata()
{
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do    
	ls -l data/run${RUNNUMBER1[$I]}_$MNSFX.v$MNINVERSION.root
    done
}


grep "CONFIG"  DoAna.sh 
anahelp

