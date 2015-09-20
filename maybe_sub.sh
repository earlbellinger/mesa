#!/bin/bash

#### Check if the queuing system is available. Otherwise, run locally. 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

if command -v condor_submit >/dev/null 2>&1
  then
    ## Create a shell script
    echo "#!/usr/bin/sh
$1
" > "condor-$1.sh"
    chmod +x "condor-$1.sh"
    
    ## Create and submit a condor job 
	echo "Universe   = vanilla
getenv     = True
Executable = condor-$1.sh
Output     = condor-$1.out
Error      = condor-$1.error
Log        = condor-$1.log

queue
" > "condor-$1.job"
    condor_submit "condor-$1.job"
    condor_wait -status -wait 3600 "condor-$fname.log"
  else
    eval "$1"
fi
