#! /bin/bash                                                                                                       

###-------------------------------------
RE1MX2=1
## 0 REAL and MIXed
## 1 REAL only
## 2 MIXed only

RUNNUMBER1=(\
"2273" "2274" "2275" "2276" "2283" "2284" \
"2285" "2286" "2288" "2289" "2291" "2310" "2311" "2314" "2315" "2320" "2322" "2323" "2324" "2325" \
"2331" "2332" "2333" "2334" "2335" "2336" "2337" "2340" "2341" "2362" "2363" "2368" "2369" "2370" \
"2371" "2372" "2373" "2374" "2375" "2378" "2379" "2380" "2381" "2382" "2383" "2384" "2385" "2386" \
"2387" "2388" "2389" "2391" "2392" "2393" "2394" "2395" "2396" "2397" "2398" "2399" "2400" "2401" \
"2402" "2429" "2432" "2433" "2434" "2437" "2438" "2439" "2440" "2442" "2453" "2461" "2462" "2463" \
"2501" "2502" "2503" "2505" "2506" "2507" "2508" "2509" \
)

#RUN=2841 VER=4.2.1 MIX=0 FLC=run2841_rf.v4.2.Psicv0 FLCS=run2841_rf.v4.5.Psis1rcv0 root flw_process4.C

VERSION=5.1.0
FLC=run2273_rf.v5.1.Psi2rtcv0   # Phi_rp
FLCS=run2273_rf.v5.1.Psis1rcv0  # Phi_sub 
     
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

