"""
JCD's original model S kernels

FORTRAN binaries in big_endian format

read(21)np,(rd(ji),ji=np,1,-1)
read(21,end=999)lll,nnn, ooo,(ak(ii),ii=np,1,-1)
"""

with open('data/ker.rho-gam1.1', 'rb') as f:
    g1k_file = f.read()

# get grid size
fmt = '>'+2*'i'
size = struct.calcsize(fmt)
N1, N2 = struct.unpack(fmt, g1k_file[:size])
g1k_file = g1k_file[size:]

# obtain mesh
fmt = '>'+N2*'f'
size = struct.calcsize(fmt)
x = np.array(struct.unpack(fmt, g1k_file[:size]))
g1k_file = g1k_file[size:]

while len(g1k_file) > 0:
    # get number of grid points... again 
    fmt = '>'+2*'i'
    size = struct.calcsize(fmt)
    out = struct.unpack(fmt, g1k_file[:size])
    g1k_file = g1k_file[size:]

    # get l and n
    fmt = '>'+2*'i'
    size = struct.calcsize(fmt)
    l, n = struct.unpack(fmt, g1k_file[:size])
    g1k_file = g1k_file[size:]

    # get o, whatever that is
    fmt = '>f'
    size = struct.calcsize(fmt)
    o, = struct.unpack(fmt, g1k_file[:size])
    g1k_file = g1k_file[size:]
    
    # get kernels
    fmt = '>'+N2*'f'
    size = struct.calcsize(fmt)
    ak = np.array(struct.unpack(fmt, g1k_file[:size]))
    g1k_file = g1k_file[size:]

