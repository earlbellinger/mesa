#!/bin/bash

change() { #param initval newval
    sed -i.bak "s/\!$1 = $2/$1 = $3/g" inlist_1.0
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

run() { # nameOfRun
    condor_submit mesa.job
    condor_wait condor.log
    cp final_profile.data "LOGS/profile_$1.data"
    cp final_profile.data.FGONG "LOGS/profile_$1.data.FGONG"
}

rerun() { # nameOfRun
    run $1
    tail -n+7 LOGS/history.data >> history.data
}

dualexp() {
    expname="M=$M;Y=$Y"
    dirname='dualexp'/$expname
    mkdir -p $dirname
    cp -r 1M_pre_ms_to_rg_template/* $dirname
    cd $dirname
    change 'initial_y' '-1' "$Y"
    change 'initial_mass' '1.0' "$M"
    
    ./mk
    run "ZAMS"
    mv LOGS/history.data .
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    rerun "solar-age"
    mv history.data LOGS
    cd -
}

for M in $(seq 0.8 0.1 1.2); do
    for Y in $(seq 0.23 0.01 0.33); do
        dualexp $M $Y &
    done
done

sleep 60
while condor_q | grep -lq bellinger; do sleep 60; done

for fgong in $(find dualexp -type f -name "*.FGONG"); do
    ./fgong2freqs.sh $fgong &
done

