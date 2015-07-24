#!/bin/bash

for exp in fgongs/*; do
    for ev_stage in "$exp"/RGB_bump; do
        for fgong in "$ev_stage"/*.FGONG; do
            ./fgong2freqs.sh $fgong &
            exit
        done
    done
done
