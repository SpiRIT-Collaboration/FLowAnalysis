# S015FlowAnalysis

This software enable us to analyize "Collective Flow" for SpiRIT-TPC.
There are 4 steps.
You can start from flw_process1.
In each step, shell script is prepared.

Installation
$ mkdir build
$ cd build
$ cmake ../ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
$ make -j
$ source config.sh

$ cd macros
$ soruce flw_process1.sh

