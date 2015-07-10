#!/bin/bash
mkdir -p exp_Y
for i in `seq 15 25`; do
#for i in `seq 15 25`; do
    mix_length=$(bc -l <<< 'scale=2; '$i'/10')
    #Y=$(bc -l <<< 'scale=2; '$i'/10')
    cp -r 1M_pre_ms_to_rg_template exp_Y/1M_$mix_length
    cd exp/1M_$mix_length
    repl="s/mixing_length_alpha = 2.1/mixing_length_alpha = $mix_length/g"
    sed -i.bak repl inlist_1.0
    ./mk
    condor_submit mesa.job
    cd ../..
done
