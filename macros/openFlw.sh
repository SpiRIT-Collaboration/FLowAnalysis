#/bin/bash

export RNF132s=2841,2843,2844,2845,2846,2848,2849,2850,2851,2852

export DBP=_rf.v4.3.0.cv0

RUN0={$RNF132s} DB0=$DBP  root openFlw.C
