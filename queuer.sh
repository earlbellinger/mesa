#!/bin/bash

change() { 
    sed -i.bak "s/\!$1 = $2/$1 = $3/g" inlist_1.0
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

run() {
    condor_submit mesa.job
    condor_wait condor.log
    cp final_profile.data "LOGS/profile_$1.data"
    cp final_profile.data.FGONG "LOGS/profile_$1.data.FGONG"
}

rerun() {
    run $1
    tail -n+7 LOGS/history.data >> history.data
}

exp() {
    dirname='exp'/$1
    mkdir -p $dirname
    parameter=$(bc -l <<< 'scale=2; '$i'/'$2)
    cp -r 1M_pre_ms_to_rg_template $dirname/$parameter
    cd $dirname/$parameter
    change "$3" "$4" "$parameter"
    
    ./mk
    run "zams"
    mv LOGS/history.data .
    change "write_profiles_flag" ".false." ".true."
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    rerun "H-exhausted"
    change "write_profiles_flag" ".true." ".false."
    change "xa_central_lower_limit(1)" "1d-3" "0"
    rerun "TP"
    
    mv history.data LOGS
    cd -
}

for i in `seq 16 26`; do exp 'alpha' '10' 'mixing_length_alpha' '2.1' &
done
for i in `seq 23 33`; do exp 'Y' '100' 'initial_y' '-1' &
done
for i in `seq 95 105`; do exp 'M' '100' 'initial_mass' '1.0' &
done
