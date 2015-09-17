import sys
import struct
import numpy as np
import os

def main(arguments):
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True,
                        help='eigenfunction file (.amde)',
                        type=str)
    parser.add_argument('-o', '--output', default='',
                        help='output directory',
                        type=str)
    parser.add_argument('-n', '--normalize', 
                        help="divide radius by max(radius)",
                        default=False, action='store_true')
    parser.add_argument('-d', '--suppress-dat', 
                        help="don't save eigenfunction files",
                        default=False, action='store_true')
    parser.add_argument('-p', '--suppress-plot', 
                        help="don't save eigenfunction plots",
                        default=False, action='store_true')
    args = parser.parse_args(arguments)
    parse_eig(eig_filename=args.input, 
              output_dir=args.output,
              normalize=args.normalize,
              save_dat=(not args.suppress_dat),
              save_plot=(not args.suppress_plot))

def parse_eig(eig_filename, output_dir='', normalize=False, 
              save_dat=True, save_plot=True):
    if (save_dat or save_plot) and output_dir != '' and \
        not os.path.exists(output_dir):
            os.makedirs(output_dir)
    
    if save_plot:
        import matplotlib as mpl 
        from matplotlib import pyplot as plt 
        mpl.rc('font', family='serif') 
        mpl.rc('text', usetex='true') 
        mpl.rc('text', dvipnghack='true') 
        mpl.rcParams.update({'font.size': 18}) 
    
    with open(eig_filename, 'rb') as f:
        eig_file = f.read()
    
    # number of mesh points and the mesh
    fmt = '<2i'
    size = struct.calcsize(fmt)
    N1, N2 = struct.unpack(fmt, eig_file[:size])
    eig_file = eig_file[size:]
    
    # read in the x-axis
    fmt = '<'+N2*'d'
    size = struct.calcsize(fmt)
    x = np.array(struct.unpack(fmt, eig_file[:size]))
    if normalize:
        x = x/max(x)
    eig_file = eig_file[size:]
    
    while len(eig_file)>4:
        # not sure what this is
        fmt = '<2i'
        size = struct.calcsize(fmt)
        out = struct.unpack(fmt, eig_file[:size])
        eig_file = eig_file[size:]
        
        # obtain mode summary
        fmt = '<'+50*'d'
        size = struct.calcsize(fmt)
        cs = struct.unpack(fmt, eig_file[:size])
        l = int(cs[17])
        n = int(cs[18])
        print("Extracting eigenfunction for mode l=%d, n=%d"%(l,n))
        eig_file = eig_file[size:]
        
        # eigenfunctions
        N = 2
        fmt = '<'+N*N2*'d'
        size = struct.calcsize(fmt)
        z = np.array(struct.unpack(fmt, eig_file[:size]))
        y1 = z[0::N]
        y2 = z[1::N]
        eigenfunction = np.vstack((x, y1, y2)).T
        eig_file = eig_file[size:]
        
        fname = 'eig_l=%d_n=%d'%(l,n)
        if save_dat:
            np.savetxt(os.path.join(output_dir, fname+'.dat'), eigenfunction)
        if save_plot:
            plt.axhspan(0, 0, ls='dashed')
            plt.plot(x, y1, 'r-', label='radial ($y_1$)')
            plt.plot(x, y2, 'b-', label='horizontal ($y_2$)')
            plt.xlabel("r/R")
            plt.ylabel("Displacement")
            plt.legend(loc='upper center')
            plt.title('Eigenfunctions for the $\ell=%d,\;n=%d$ mode'%(l,n))
            plt.ylim([min(0, min(y1), min(y2))-0.1,
                      max(1, max(y1), max(y2))+0.1])
            if normalize:
                plt.xlim([0, 1.01])
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, fname+'.png'))
            plt.close()
    
    if len(eig_file) != 4:
        print("Error: failed to parse eigenfunction")

if __name__ == '__main__':
    exit(main(sys.argv[1:]))

