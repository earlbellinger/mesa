#!/bin/bash

## Run mix length experiments
dirname=exp_alpha
mkdir -p $dirname
for i in `seq 16 26`; do
    mix_length=$(bc -l <<< 'scale=2; '$i'/10')
    cp -r 1M_pre_ms_to_rg_template $dirname/1M_$mix_length
    cd $dirname/1M_$mix_length
    repl="s/mixing_length_alpha = 2.1/mixing_length_alpha = $mix_length/g"
    sed -i.bak "$repl" inlist_1.0
    ./mk
    condor_submit mesa.job
    cd ../..
done

## Run helium experiments
dirname=exp_Y
mkdir -p $dirname
#for i in `seq 230 5 300`; do
for i in `seq 23 33`; do
    #Y=$(bc -l <<< 'scale=3; '$i'/1000')
    Y=$(bc -l <<< 'scale=2; '$i'/100')
    cp -r 1M_pre_ms_to_rg_template $dirname/1M_$Y
    cd $dirname/1M_$Y
    repl="s/initial_y = -1/initial_y = $Y/g"
    sed -i.bak "$repl" inlist_1.0
    ./mk
    condor_submit mesa.job
    cd ../..
done
