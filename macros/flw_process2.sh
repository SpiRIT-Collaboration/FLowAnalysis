#! /bin/bash                 

######>>>  RUN=2841 VER=4.1 MIX=0 root flw_process2.C\(2000\) 

source runList.sh

source setup.sh
##--- for Process2 ------------------------------------
# *****> <Edit Here>  

#RUNNUMBER1=("2273")    
#RUNNUMBER1=(${RNF132})    
#RUNNUMBER1=(${RNF108})                                                                              
#RUNNUMBER1=(${RNF124})
#RUNNUMBER1=(${RNF112})
RUNNUMBER1=(${RNF132r} ${RNF108}) # ${RNF124} ${RNF112})    
# *****> <Edit Here>                                                                                                            
VERSION=12.0

# *****> <Edit Here> 
RE1MX2=1
## 0 REAL and MIXed
## 1 REAL only
## 2 MIXed only

function process2() {
    typeset -i I=0
    
    if [ ${RE1MX2} -ne "2" ]; then
	MIX=0
	echo " REAL data --- "

	while(( I < ${#RUNNUMBER1[@]} ))
	do
            RUN=${RUNNUMBER1[I]}
	    LOG=log/proc2_${RUN}_v${VERSION}_${MIX}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} root flw_process2.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} root -b -q flw_process2.C >& $LOG &
            let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi
	    
            RUN=${RUNNUMBER1[I]}
	    LOG=log/proc2_${RUN}_v${VERSION}_${MIX}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} root flw_process2.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} root -b -q flw_process2.C >& $LOG &
            let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi
	    
            RUN=${RUNNUMBER1[I]}
	    LOG=log/proc2_${RUN}_v${VERSION}_${MIX}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} root flw_process2.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} root -b -q flw_process2.C >& $LOG 
            let I++
	done
    fi

    if [ ${RE1MX2} -ne "1" ]; then
	echo " MIXed data --- "
	MIX=1
	I=0
	while(( I < ${#RUNNUMBER1[@]} ))
	do
	    RUN=${RUNNUMBER1[I]}
	    LOG=log/proc2_${RUN}_v${VERSION}_${MIX}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} root flw_process2.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} root -b -q flw_process2.C >& $LOG &
	    let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi

	    RUN=${RUNNUMBER1[I]}
	    LOG=log/proc2_${RUN}_v${VERSION}_${MIX}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} root flw_process2.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} root -b -q flw_process2.C >& $LOG &
	    let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi
	    
	    RUN=${RUNNUMBER1[I]}
	    LOG=log/proc2_${RUN}_v${VERSION}_${MIX}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} root flw_process2.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} root -b -q flw_process2.C >& $LOG 
	    let I++
	done
    fi

}


process2

