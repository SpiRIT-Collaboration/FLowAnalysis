#! /bin/bash
#
# mc_process1.sh
# (c) Mizuki Nishimura
#
# ----------------------

# MC data

export MCDIR=/cache/scr/spirit/kaneko/rootfiles/recoDataMC/20180711dev_20180713dev.RC/
export STVERSION=1660.2c1f3c3
export MCSYSTEM=amd_132Sn124Sn270AMeV_cluster_SLy4

##--- for Process1 ------------------------------------
# *****> <Edit Here>

# *****> <Edit Here>     
VERSION=0

source setup.sh


echo RUN=${RUN} VER=$VERSION root -b -q flw_process1.C '>&' $LOG
RUN=${RUN} VER=$VERSION root -b -q flw_process1.C >& $LOG 







