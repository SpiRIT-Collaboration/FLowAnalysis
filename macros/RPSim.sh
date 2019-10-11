#!/bin/bash

source setup.sh
root -b -q RPSim.C\(10,100000\) &
root -b -q RPSim.C\(11,100000\) &
root -b -q RPSim.C\(12,100000\) 
