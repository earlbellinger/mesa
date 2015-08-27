#!/bin/bash

for exp in fgongs/*; do
    #for ev_stage in $exp/RGB; do
    for ev_stage in $exp/*; do
        for fgong in "$ev_stage"/*.FGONG; do
            if grep -q 'RGB' <<< $fgong; then continue; fi
            if [[ -f $fgong ]]; then
                ./fgong2freqs.sh $fgong &
                #echo $fgong
            fi
        done
    done
done

