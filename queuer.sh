#!/bin/bash

change() {
    sed -i.bak "s/\!$1 = $2/$1 = $3/g" inlist_1.0
    sed -i.bak "s/$1 = $2/$1 = $3/g" inlist_1.0
}

run() {
    condor_submit mesa.job
    condor_wait condor.log
    mv final_profile.data "LOGS/profile_$1.data"
    mv final_profile.data.FGONG "LOGS/profile_$1.data.FGONG"
}

rerun() {
    run $1
    tail -n+7 LOGS/history.data >> history.data
}

exp() { # dirname scale divby inlistname defaultval
    dirname=$1
    paramname=$2
    initval=$3
    newval=$4
    mkdir -p $dirname
    if [ ${newval:0:1} == "." ]; then newval='0'$newval; fi
    cp -r 1M_pre_ms_to_rg_template $dirname/$newval
    cd $dirname/$newval
    change "$paramname" "$initval" "$newval"
    
    ./mk
    run "ZAMS"
    mv LOGS/history.data .
    ##change "write_profiles_flag" ".false." ".true."
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    rerun "solar-age"
    #change "max_age" "4.57e9" "-1"
    #rerun "H-exhausted"
    ##change "xa_central_lower_limit(1)" "1d-3" "0"
    ##rerun "RGB"
    
    mv history.data LOGS
    cd -
}

for i in $(seq 1.6 0.025 2.6)
  do exp 'exp/alpha' 'mixing_length_alpha' '2.1' $i &
done

for i in $(seq 0.9 0.0025 1.1)
  do exp 'exp/M' 'initial_mass' '1.0' $i &
done

for i in $(seq 0.23 0.00125 0.33)
  do exp 'exp/Y' 'initial_y' '-1' $i &
done

#for i in `seq 44 124`; do exp 'alpha' '3' '40' 'mixing_length_alpha' '2.1' &
#done
#for i in `seq 184 264`; do exp 'Y' '5' '800' 'initial_y' '-1' &
#done
#for i in `seq 360 440`; do exp 'M' '4' '400' 'initial_mass' '1.0' &
#done
