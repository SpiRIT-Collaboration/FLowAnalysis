#! /Bin/Bash
#
# DoAna.Sh
# (C) Mizuki Nishimura
#
# ----------------------

source setup.sh
source runList.sh
export MNMACRO=DoAna_phys.C          ##<--- Flow analysis MACRO name

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

export MDCUT=0. #0.0               ##   <------@@ mid-rapidity cut
##------>>> EDIT HERE 
#export MNLC=0
#export MNUC=50

multL=55,46,0,46,28,30,46
multU=80,55,46,80,49,48,54
##     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17
#####" 0,20,35,40,50,55,42,42,40,42, 0, 0,52,47,46,51,52,42,47"
#####"35,40,45,55,65,80,52,56,54,55,54,55,55,55,54,80,80,80,80"
###export MNTRK=13
###export MNTRK=15

export MNTRK=0    #Central
export PHICUT=2   #|phi|>140
export MNSUBVERSION=15
#----///
export MNTRK=5    #Tommy mid-central for 132
export PHICUT=6   #-30<phi<20or|phi|>140 
export MNSUBVERSION=19
#----///
#----///
export MNTRK=4    #Tommy mid-central for 108
export PHICUT=6   #-30<phi<20or|phi|>140 
export MNSUBVERSION=19
#----///
#--/
export MNTRK=0    #central
export PHICUT=5   #-30<phi<20
export MNSUBVERSION=18
#----///
#--/
export MNTRK=6    #Tommy2 mid-central
export PHICUT=5   #-30<phi<20or|phi|>140 
export MNSUBVERSION=23
#--/
#--/
export MNTRK=6    #Tommy2 mid-central
export PHICUT=5   #-30<phi<20or|phi|>140 
export MNSUBVERSION=24
#--/
#--/
export MNTRK=1    #mid-central
export PHICUT=5   #-30<phi<20
export MNSUBVERSION=0
#----///
#--/
export MNTRK=1    #mid-central
export PHICUT=6   #-30<phi<20or|phi|>140 
export MNSUBVERSION=1
#--/
export MNTRK=1    #mid-central
export PHICUT=5   #-30<phi<20
export MNSUBVERSION=4
#----///
#--/
export MNTRK=1    #mid-central
export PHICUT=5   #-30<phi<20
export MNSUBVERSION=6
#----///
#--/
export MNTRK=1    #mid-central
export PHICUT=5   #-30<phi<20
export MNSUBVERSION=5
#----///
#--/
export MNTRK=1    #mid-central
export PHICUT=5   #-30<phi<20
export MNSUBVERSION=0
#----///
export MNTRK=1    #mid-central
export PHICUT=6   #-30<phi<20or|phi|>140  
export MNSUBVERSION=1
#--/
#--/
export MNTRK=6    #Tommy2 mid-central
export PHICUT=6   #-30<phi<20or|phi|>140 
export MNSUBVERSION=2
#--/

#########----------------------
MLC=(${multL})
MUC=(${multU})
export MNLC=${MLC[$MNTRK]}
export MNUC=${MUC[$MNTRK]}

export MNSFX=flow
export MNINVERSION=3             ##   <------@@ input
export MNOUTVERSION=3.0            ##   <------@@ output 

function doflow3fmFull()
{
    doflow3fm112
    doflow3fm124
    doflow3fm108
    doflow3fm132
}

function doflow3fm108s() 
{
    data108s
    PARTICLES=("0")
    doflowmulti $MNSUBVERSION
}
function doflow3fm112one() 
{
    data112one
    doflow 2 $MNSUBVERSION
}
function doflow3fm112p() 
{
    data112
    doflow 2 $MNSUBVERSION
}
function doflow3fm132() 
{
    data132
    PARTICLES=("0" "1" "2" "3" "4")
#    PARTICLES=("0")
    doflowmulti $MNSUBVERSION
}
function doflow3fm108() 
{
    data108
#    PARTICLES=("0" "1" "2" "3" "4")
    PARTICLES=("0")
    doflowmulti $MNSUBVERSION
}
function doflow3fm112() 
{
    data112
    PARTICLES=("0" "1" "2" "3" "4" "5" "6")
#    PARTICLES=("5" "6")
    doflowmulti $MNSUBVERSION
}
function doflow3fm124() 
{
    data124
    PARTICLES=("0" "1" "2" "3" "4" "5" "6")
    doflowmulti $MNSUBVERSION
}
##--
#
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

function data108one()    ##CONFIG
{
##--
    export MNSN=108
    export MNRNF=$RNF108one
    export MNOSUBV=

    if [ -n "$1" ]; then
	export MNINVERSION=$1            ##   <------@@ input
	export MNOUTVERSION=$2           ##   <------@@ output 
    fi

    commonsetup
##<----
}

function data112one()    ##CONFIG
{
##--
    export MNSN=112
    export MNRNF=$RNF112one
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
    if [ $MNTRK -eq 10 ]; then
	export MNTRK=14
    fi
    if [ $MNTRK -eq 11 ]; then
	export MNTRK=15
    fi

    export MNLC=${MLC[$MNTRK]}
    export MNUC=${MUC[$MNTRK]}
    export MNSN=132
    export MNRNF=$RNF132
    export MNRNF=$RNF132_BH
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


function doflowmultiHeavy()
{
    PARTICLES=("3" "4" "5" "6")
    doflowmulti $1
}

function doflowmultiall()
{
    PARTICLES=("0" "1" "2" "3" "4" "5" "6")
    doflowmulti $1
}

function doflowmulti()
{
    if [ -n "$1" ]; then
	export OSUBV=$1
	echo "output version -> "  $OSUBV
    fi


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
    echo " Step 3 > doflow       :: open files "
    echo "          doflow 2  0(output version#) 1(Corretion versoin)::" $MNMACRO
    echo "                  -1:pi- 1:pi+, 2:p, 3:d, 4:t, 5:3He, 6: 4He, 7:n, 8:H" 
    echo "                  run #(partid) #(Output version)"
    echo "++++++++++++++++++++++++++++++++++++++++"
    echo "Change runnumber : export MNRNF=\$RNF132s then commonsetup"
    echo "------ Setup for temporarly ------------"
    echo export MNSFX=$MNSFX
    echo export MNINVERSION=$MNINVERSION    ##   <------@@ input
    echo export MNOUTVERSION=$MNOUTVERSION  ##   <------@@ output 
    echo "Then >commonsetup "
    echo "----------------------------------------"
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
    echo "Output subVer    ::"  $MNSUBVERSION
    echo "Selected track   ::"  $TRK
    echo "Number of event  ::"  $MXEVT
    echo "Mid-Y cut        ::"  $MDCUT
    echo "Min multiplicity ::"  $LC
    echo "Max maltiplicity ::"  $UC
    echo "Phi Angle cut    ::"  $PHICUT

    echo "++++++++++++++++++++++++++++++++++++++++"
    echo export MNINVERSION=$IVER            ##   <------@@ input
    echo export MNOUTVERSION=$OVER           ##   <------@@ output 
    echo "++++++++++++++++++++++++++++++++++++++++"

}

function checkdata()
{
    typeset -i I=0
    while(( $I < ${#RUNNUMBER1[@]} ))
    do    
	ls -l data/run${RUNNUMBER1[I]}_$MNSFX.v$MNINVERSION.root
#	echo $I	
#	ls -l data/run${RUNNUMBER1[I]}_'s_Proton'.v$MNINVERSION.root
	let I++
    done
}

#showenv
#grep "function" DoAna.sh
grep "CONFIG"   DoAna.sh 
anahelp
