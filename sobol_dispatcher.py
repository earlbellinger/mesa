#!usr/bin/env python

from sobol_lib import i4_sobol
import numpy as np
import subprocess
from time import sleep
from math import log10

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-M', default=[0.8, 1.2], nargs=2, type=float,
                        help='range of masses')
    parser.add_argument('-Y', default=[0.22, 0.32], nargs=2, type=float, 
                        help='range of helium values')
    parser.add_argument('-Z', default=[0.0001, 0.04], nargs=2, type=float,
                        help='range of metallicity values')
    parser.add_argument('-a', '--alpha', default=[1.7, 2.2], nargs=2,type=float,
                        help='range of mixing length parameter values')
    parser.add_argument('-N', default=100, help='number of models to generate',
                        type=int)
    parser.add_argument('-s', '--seed', default=0, type=int,
                        help='offset for sobol numbers')
    args = parser.parse_args(arguments)
    print(args)
    dispatch(ranges=np.vstack((args.M, args.Y, 10**np.array(args.Z), args.alpha)),
             N=args.N, logs=[0, 0, 1, 0], seed=args.seed)

def dispatch(ranges, N, logs, seed=0):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    for i in range(seed, N+seed):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        for j, val in enumerate(vals):
            if (logs[j]):
                vals[j] = log10(val)
        bash_cmd = "maybe_sub.sh dispatch.sh -M %.6f -Y %.6f -Z %.6f -a %.6f"%\
            tuple(val for val in vals)
        subprocess.Popen(bash_cmd.split(), shell=False)
        sleep(0.2)

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))
