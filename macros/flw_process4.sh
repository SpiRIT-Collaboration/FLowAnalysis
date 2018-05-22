#! /bin/bash                                                                                                       

###-------------------------------------
RE1MX2=1
## 0 REAL and MIXed
## 1 REAL only
## 2 MIXed only

RUNNUMBER1=(\
"2542" "2543" "2544" "2545" "2546" "2547" "2548" "2549" \
"2550" "2551" "2552" "2553" "2554" "2555" "2556" "2557" "2558" "2559" \
"2560" "2562" "2563" "2564" "2565" "2566" "2567" "2568" "2569" "2570" \
"2571" "2572" "2573" "2574" "2575" "2576" "2578" "2579" "2580" "2581" \
"2582" "2583" "2584" "2585" "2586" "2587" "2588" "2589" "2590" "2591" \
"2592" "2593" "2594" "2595" "2596" "2597" "2598" "2599" "2600" "2601" \
"2603" "2604" "2605" "2607" "2609" "2610" \
"2613" "2617" "2618" "2619" "2620" "2621" "2622" "2623" \
"2624" "2632" "2633" "2634" "2643" "2645" "2646" "2647" "2648" \
"2649" "2650" "2652" "2653" \
)

#RUN=2841 VER=4.2.1 MIX=0 FLC=run2841_rf.v4.2.Psicv0 FLCS=run2841_rf.v4.5.Psis1rcv0 root flw_process4.C

VERSION=5.0.0
FLC=run2542_rf.v5.0.Psi2rtcv0   # Phi_rp
FLCS=run2542_rf.v5.0.Psis1rcv0  # Phi_sub 
     
echo $RUNNUMBER1

function process4() {
    typeset -i I=0
    
    if [ ${RE1MX2} -ne "2" ]; then
	MIX=0
	echo " REAL data --- "

	while(( I < ${#RUNNUMBER1[@]} ))
	do
            RUN=${RUNNUMBER1[I]}
	    LOG=log/prooc4_${RUN}_v${VERSION}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root flw_process4.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root -b -q flw_process4.C >& $LOG &
            let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi
	    
            RUN=${RUNNUMBER1[I]}
	    LOG=log/prooc4_${RUN}_v${VERSION}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root flw_process4.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root -b -q flw_process4.C >& $LOG &
            let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi
	    
            RUN=${RUNNUMBER1[I]}
	    LOG=log/prooc4_${RUN}_v${VERSION}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root flw_process4.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root -b -q flw_process4.C >& $LOG 
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
	    LOG=log/prooc4_${RUN}_v${VERSION}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root flw_process4.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root -b -q flw_process4.C >& $LOG &
	    let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi

	    RUN=${RUNNUMBER1[I]}
	    LOG=log/prooc4_${RUN}_v${VERSION}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root flw_process4.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root -b -q flw_process4.C >& $LOG &
	    let I++
	    if [ $I -ge ${#RUNNUMBER1[@]} ]; then
		break;
	    fi
	    
	    RUN=${RUNNUMBER1[I]}
	    LOG=log/prooc4_${RUN}_v${VERSION}.log
	    echo RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root flw_process4.C  '>&' $LOG 
	    RUN=${RUN} VER=${VERSION} MIX=${MIX} FLC=${FLC} FLCS=${FLCS} root -b -q flw_process4.C >& $LOG 
	    let I++
	done
    fi

}


process4

