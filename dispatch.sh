#!/bin/bash

dualexp() { 
    expname="M=$M""_""Y=$Y""_""Z=$Z""_""alpha=$alpha"
    dirname="dexp/$expname"
    mkdir -p "$dirname"
    cp -r dualexp_template/* "$dirname"
    cd "$dirname"
    change 'initial_mass' '1.0' "$M"
    change 'initial_y' '-1' "$Y"
    change 'initial_z' '0.02' "$Z"
    change 'Zbase' '0.02' "$Z"
    change 'mixing_length_alpha' '2.1' "$alpha"
    
    ./mk
    run "ZAMS"
    mv LOGS/history.data .
    change "stop_near_zams" ".true." ".false."
    change "create_pre_main_sequence_model" ".true." ".false."
    change "load_saved_model" ".false." ".true."
    change "write_profiles_flag" ".false." ".true."
    rerun "Hexh"
    mv history.data LOGS
    
    cd -
    
    for fgong in $(find "$dirname/LOGS" -type f -name "*.FGONG"); do
        getfreqs "$fgong" &
        sleep 1
    done
}

getfreqs() { # fgongfile
    for i in 1 2; do 
        ./fgong2freqs.sh "$1" && break || sleep 60
    done
}

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

while [ "$#" -gt 0 ]; do
  case "$1" in
    -M) M="$2"; shift 2;;
    -Y) Y="$2"; shift 2;;
    -Z) Z="$2"; shift 2;;
    -a) alpha="$2"; shift 2;;

    --mass=*) M="${1#*=}"; shift 1;;
    --helium=*) Y="${1#*=}"; shift 1;;
    --metal=*) Z="${1#*=}"; shift 1;;
    --alpha=*) alpha="${1#*=}"; shift 1;;
    --mass|--helium|--metal|--alpha) echo "$1 needs argument" >&2; exit 1;;

    -*) echo "unknown option: $1" >&2; exit 1;;
    *) handle_argument "$1"; shift 1;;
  esac
done

dualexp


