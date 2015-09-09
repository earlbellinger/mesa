#!/bin/bash

#### Converter for FGONG to oscillation mode frequencies using ADIPLS 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

### Parse command line tokens 
## -h and --help will display help instructions 
if [ -z ${1+x} ] || [ $1 == '-h' ] || [ $1 == '--help' ]; then
    echo "Converter for FGONG to oscillation mode frequencies using ADIPLS."
    echo "Usage: ./fgong2freqs.sh FGONG OUTPUT"
    echo
    echo "Will make a directory called OUTPUT and place OUTPUT.dat in there."
    echo "In absence of OUTPUT, will be created from the FGONG filename."
    exit
fi

## Check that the first input (FGONG file) exists
if [ ! -e $1 ]; then
    echo "Error: Cannot locate FGONG file" $1
    exit 1
fi

## Pull out the name of the FGONG file
fname="$(basename $1)"
fname="${fname%%.*}-freqs"

## If the second (OUTPUT) argument doesn't exist, create one from the first
if [ -z ${2+x} ]; then
    path=$(dirname $1)/$fname
  else
    path=$2
fi

## Create a directory for the results and go there
## Also convert the FGONG file to AMDL format and put it in the new directory
mkdir "$path"
#fname=(${path//\/}) # remove slashes 
fgong-amdl.d $1 "$path/$fname.amdl"
cd "$path"

logfile="fgong2freqs.log"
exec > $logfile 2>&1

## Check that the amdl file was created 
if [ ! -e "$fname.amdl" ]
  then
    echo "Error: Conversion of $1 to .amdl format failed"
    exit 1
fi

## Create a redistribute file and pass it to ADIPLS' redistrb 
echo "
2 '$fname.amdl'    @    
3 '$fname-6202'    @
-1 ''        @
nn,icnmsh
6202,,,  @
icase,icvzbn,nsmth,ndisc,dlxdsc,dlgrmx,cacvzb
21   ,      ,     ,200  ,      ,5.    ,      ,@  
cg,  cx,ca  ,sig1,sig2,lmax,alphsf,adda4,accrm
1.,0.05,0.05,0.  ,0.  ,   2,      ,0.02 ,0.01,  @
nout,cn,irsu,unew   
60,,,,,,,,,  @
nmodel,kmodel,itsaml,ioldex
,,,,,,,,,,,  @
" > "redistrb-$fname.in"
redistrb.c.d "redistrb-$fname.in"

## Check that the redistribution was successful
if [ ! -e "$fname-6202" ]; then
    echo "Error: Redistribution of $fname failed"
    exit
fi

## Create an adipls.in file with some decent (?) settings 
echo "
 2  '$fname-6202'   @
 9  '$fname.log'   @
 11 '$fname.agsm'   @
 4  '$fname.amde'   @
 15 '$fname.ssm'    @
 -1 ''   @
  cntrd,
mod.osc.cst.int.out.dgn     @

mod:
  ifind,xmod,imlds,in,irname,nprmod,
       ,    ,     ,  ,      ,      ,  @
  ntrnct,ntrnsf,imdmod,
       ,       ,      , @
osc:
  el,nsel,els1,dels,dfsig1,dfsig2,nsig1,nsig2,
    ,    4,  0,   1,    0,       ,    1,    1,   @
  itrsig,sig1,istsig,inomde,itrds,
       1,   5,     1,     1,    , @
  dfsig,nsig,iscan,sig2,
       ,  2,   5000, 10000,     @
  eltrw1, eltrw2, sgtrw1, sgtrw2,
       0,     -1,     0,       -1,  @
cst:
  cgrav
  6.672320e-8               @
int:
  iplneq,iturpr,icow,alb,
  0,0,0,1,,,,,,             @
  istsbc,fctsbc,ibotbc,fcttbc,
  1,0,0,0,,,,,  @
  mdintg,iriche,xfit,fcnorm,eps,epssol,itmax,dsigre,
  5,1,0.99,,,,15,0,,,,,,  @
  fsig,dsigmx,irsevn,xmnevn,nftmax,itsord
  0.001,0.1,2,0,1,20,,,,,  @
out:
  istdpr,nout,nprcen,irsord,iekinr
  9,50,100,20,0,,,,,,,,     @
  iper,ivarf,kvarf,npvarf,nfmode,
  1,1,2,0,3,,,,,,,,     @
  irotkr,nprtkr,igm1kr,npgmkr,ispcpr,
  0,,0,,0,,,,,,,     @
  icaswn, sigwn1, sigwn2, frqwn1, frqwn2,iorwn1, iorwn2, frlwn1, frlwn2
        ,      0,     -1,      0,  10000,  -200,     50,      0,     -1,   @
dgn:
  itssol,idgtss,moddet,iprdet,npout
  0,0,0,0,,,,,,,,,,,     @
  imstsl,imissl,imjssl,idgnrk
  ,,,,,,,,,    @
" > "adipls-$fname.in"

### Time to Run ADIPLS!
## Check if the queuing system is available. Otherwise, run locally.
if command -v condor_submit >/dev/null 2>&1
  then
    ## Create a shell script for executing adipls
    echo "#!/usr/bin/sh
adipls.c.d adipls-$fname.in
" > "adipls-$fname.sh"
    chmod +x "adipls-$fname.sh"
    
    ## Create and submit a condor job for executing that shell script
	echo "Universe   = vanilla
getenv     = True
Executable = adipls-$fname.sh
Output     = condor-$fname.out
Error      = condor-$fname.error
Log        = condor-$fname.log

queue
" > "condor-$fname.job"
    condor_submit "condor-$fname.job"
    condor_wait -status -wait 3600 "condor-$fname.log"
  else
    adipls.c.d "adipls-$fname.in"
fi

## Check that the .agsm file was created successfully, and convert it to freqs 
if [ ! -e "$fname.agsm" ]; then
    echo "Error: Failed to generate $fname.agsm"
    exit 1
fi
set-obs.d 6 "$fname.agsm" "$fname.dat" 3

## Check that the frequencies were created, and remove the annoying text 
if [ ! -e "$fname.dat" ]; then
    echo "Error: Failed to generate frequency information $fname.dat"
    exit 1
fi
cp "$fname.dat" "$fname.dat.bak"
cat "$fname.dat.bak" | tr -s ' ' | cut -d ' ' -f 2-4 > "$fname.dat"
#sed -i.bak 's/  [a-zA-Z].*$//g' "$fname.dat"

### Hooray!
cp "$fname.dat" ..
echo "Conversion complete. Results can be found in $path/$fname.dat"
exit 0

