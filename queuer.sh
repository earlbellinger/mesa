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

exp() { # parameter scale divby inlistname defaultval
    dirname='exp'/$1
    mkdir -p $dirname
    parameter=$(bc -l <<< 'scale='$2'; '$i'/'$3)
    if [ ${parameter:0:1} == "." ]; then parameter='0'$parameter; fi
    cp -r 1M_pre_ms_to_rg_template $dirname/$parameter
    cd $dirname/$parameter
    change "$4" "$5" "$parameter"
    
    ./mk
    run "ZAMS"
    mv LOGS/history.data .
    #change "write_profiles_flag" ".false." ".true."
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    rerun "solar-age"
    change "max_age" "4.57e9" "-1"
    rerun "H-exhausted"
    #change "xa_central_lower_limit(1)" "1d-3" "0"
    #rerun "RGB"
    
    mv history.data LOGS
    cd -
}

for i in `seq 44 124`; do exp 'alpha' '3' '40' 'mixing_length_alpha' '2.1' &
done
for i in `seq 184 264`; do exp 'Y' '5' '800' 'initial_y' '-1' &
done
for i in `seq 360 440`; do exp 'M' '4' '400' 'initial_mass' '1.0' &
done

