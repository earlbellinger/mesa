#!usr/bin/env python

from sobol_lib import i4_sobol
import numpy as np
import subprocess
from math import log10
from time import sleep

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-M', default=[0.8, 2], nargs=2,
                        help='range of masses')
    parser.add_argument('-Y', default=[0.2, 0.4], nargs=2,
                        help='range of helium values')
    parser.add_argument('-Z', default=[1.002, 1.071], nargs=2,
                        help='logarithmic range of helium values')
    parser.add_argument('-a', '--alpha', default=[1.5, 2.5], nargs=2,
                        help='range of mixing length parameter values')
    parser.add_argument('-N', default=100, help='number of models to generate',
                        type=int)
    parser.add_argument('-s', '--seed', default=0, type=int,
                        help='offset for sobol numbers')
    args = parser.parse_args(arguments)
    dispatch(ranges=np.vstack((args.M, args.Y, args.Z, args.alpha)), 
             logs=[0, 0, 1, 0],
             N=args.N, seed=args.seed)

def dispatch(ranges, logs, N, seed=0):
    shift = ranges[:,0]
    scale = np.array([(b-a) for a,b in ranges])
    for i in range(seed, N+seed):
        vals = shift+np.array(i4_sobol(len(ranges), i)[0])*scale
        for j, val in enumerate(vals):
            if (logs[j]):
                vals[j] = log10(val)
        bashCommand = "dispatch.sh -M %.5g -Y %.5g -Z %.5g -a %.5g" % \
            tuple(val for val in vals)
        #print(bashCommand.split())
        subprocess.Popen(bashCommand.split(), shell=False)
        sleep(1)

if __name__ == '__main__':
    import sys
    exit(main(sys.argv[1:]))

