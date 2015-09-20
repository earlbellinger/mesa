#!/bin/bash
for i in $(seq 0 50 450); do 
    python3 sobol_dispatcher.py -N $(($i+50)) -s $i
    sleep 3600
done
