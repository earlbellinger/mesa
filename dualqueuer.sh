#!/bin/bash

change() { #param initval newval
    sed -i.bak "s/\!$1 = $2/$1 = $3/g" inlist_1.0
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

run() { # nameOfRun
    condor_submit mesa.job
    condor_wait condor.log
    mv final_profile.data "LOGS/profile_$1.data"
    mv final_profile.data.FGONG "LOGS/profile_$1.data.FGONG"
}

rerun() { # nameOfRun
    run $1
    tail -n+7 LOGS/history.data >> history.data
}

dualexp() { # mass helium metallicity
    expname="M=$1;Y=$2;Z=$3"
    dirname='dualexp'/$expname
    mkdir -p $dirname
    cp -r dualexp_template/* $dirname
    cd $dirname
    change 'initial_mass' '1.0' "$1"
    change 'initial_y' '-1' "$2"
    change 'initial_z' '0.02' "$3"
    change 'Zbase' '0.02' "$3"
    
    ./mk
    run "ZAMS"
    mv LOGS/history.data .
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    change "write_profiles_flag" ".false." ".true."
    rerun "H-exhausted"
    mv history.data LOGS
    cd -
}

for M in $(seq 0.8 0.1 1.2); do
    for Y in $(seq 0.23 0.01 0.33); do
        for Z in 0.001 0.004 0.01 0.02 0.03; do
            dualexp $M $Y $Z &
        done
    done
done

sleep 60
while condor_q | grep -lq bellinger; do sleep 60; done

for fgong in $(find dualexp -type f -name "*.FGONG"); do
    ./fgong2freqs.sh $fgong &
done

