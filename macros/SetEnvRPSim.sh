#! /bin/bash

if [ $MNVERSION -eq 3 ]; then
    export MNRNF="0014,0015,0016,0020,0021"  
    export MNDBVERSION=2
elif [ $MNVERSION -eq 5 ]; then
    export MNRNF="0100"       
elif [ $MNVERSION -eq 6 ]; then
    export MNRNF="0019"       
elif [ $MNVERSION -eq 7 ]; then
    export MNRNF="0110"       
elif [ $MNVERSION -eq 8 ]; then
    export MNRNF="0110"       
elif [ $MNVERSION -eq 9 ]; then
    export MNRNF="0023"       
elif [ $MNVERSION -eq 10 ]; then
    export MNRNF="0024"       
elif [ $MNVERSION -eq 11 ]; then
    export MNRNF="0026"       
elif [ $MNVERSION -eq 12 ]; then
    export MNRNF="0025"       
elif [ $MNVERSION -eq 13 ]; then
    export MNRNF="0027"       
elif [ $MNVERSION -eq 14 ]; then
    export MNRNF="0111"       
elif [ $MNVERSION -eq 15 ]; then
    export MNRNF="0112"       
elif [ $MNVERSION -eq 16 ]; then
    export MNRNF="0028"       
elif [ $MNVERSION -eq 17 ]; then
    export MNRNF="0029"       
elif [ $MNVERSION -eq 18 ]; then
    export MNRNF="0030"       
elif [ $MNVERSION -eq 19 ]; then
    export MNRNF="0031"       
elif [ $MNVERSION -eq 20 ]; then
    export MNRNF="0032"       
elif [ $MNVERSION -eq 21 ]; then
    export MNRNF="0033"       
fi

echo run=$MNRNF