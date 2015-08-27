#!/bin/bash

Rscript 16Cyg.R &

Rscript modelS.R &

./queuer.sh

sleep 600
while condor_q | grep -lq bellinger; do sleep 60; done

Rscript HR.R

./adiplser.sh

sleep 600
while condor_q | grep -lq bellinger; do sleep 60; done

Rscript fgong.R 

