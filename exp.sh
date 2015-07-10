#!/bin/bash

mkdir exp
for fgong in ../MESA/exp/*/LOGS/profile275.data.FGONG; do
    output=$(echo $fgong | sed -e "s/..\/MESA\/exp\/1M_//g" \
                         | sed -e "s/\/LOGS\/profile275.data.FGONG//g")
    ./fgong2freqs.sh "$fgong" "exp/$output" &
done
