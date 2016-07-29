import os
import random
import sys
import time


# Converts from index to 3D cartesian coords.
def ind2sub(index, dims):
    z = index % dims[2]
    y = (index / dims[2]) % dims[1]
    x = (index / (dims[2]*dims[1])) % dims[0]
    return int(x), int(y), int(z)


# Takes a subvolume (x, y, z) index and returns
# a random position in space between (x, y, z)
# and (x+1, y+1, z+1).
def fuzzpos(sbv):
    return [i + random.random() for i in sbv]


# Writes to the given file the vtk header for a
# dataset containing npoints molecules.
def writeheader(f, npoints):
    f.write("# vtk DataFile Version 2.0\n")
    f.write("Example cartesian space datafile\n")
    f.write("ASCII\nDATASET POLYDATA\n")
    f.write("POINTS " + str(npoints) + " FLOAT\n")

    
# Writes to the given file all the points of
# specie spi.
def writepoints(f, subv, spi):
    for i in subv:
        for s in range(subv[i][spi]):
            f.write("  ")
            x, y, z = fuzzpos((i[0], i[1], i[2]))
            f.write(str(x) + " ")
            f.write(str(y) + " ")
            f.write(str(z) + "\n")    


# Writes to the given file the dataset that
# colors the molecules of the different species.
def writecolors(f, spnum):
    f.write("SCALARS Specie_" + str(sum(spnum)) + " float\n")
    f.write("LOOKUP_TABLE default\n")
    for i in range(len(spnum)): # for all specie
        for c in range(spnum[i]): # for every point in that specie
            f.write("  " + str(i) + ".0\n") # write a dataset line


            
if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: state2vtk.py <statefile> xdim ydim zdim")
        sys.exit(1)

    print("\n  --------------------------------")
    print("  |   system-state .vtk writer   |")
    print("  --------------------------------\n")
        
    xdim, ydim, zdim = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

    # open and parse .dat system state file
    with open(sys.argv[1]) as f:
        insize = os.path.getsize(sys.argv[1])

        print("    parsing .dat file.. ".format(insize/1024.0), \
            end="")

        startt = time.clock()
        
        # read 'subvolumes=' line
        t, sbc = f.readline().strip('\n').split('=')
        if t != "subvolumes":
            print("line 1: no 'subvolumes' keyword")
            sys.exit(1)
        sbc = int(sbc)
        
        # read 'species=' line
        t, spc = f.readline().strip('\n').split('=')
        if t != "species":
            print("line 2: no 'species' keyword")
            sys.exit(1)
        spc = int(spc)
        
        # read 'time=' line
        t, timestamp = f.readline().strip('\n').split('=')
        if t != "time":
            print("line 3: no 'time' keyword")
            sys.exit(1)
        timestamp = float(timestamp)

        f.readline() # eat newline

        # allocate subvolumes dict
        subv = {i:[] for i in range(sbc)}

        # read subvolumes lines
        for i in range(sbc):
            subv[i] = [int(n) for n in f.readline().split()]

        # convert points from index to subscript representation
        subv = {ind2sub(i, [xdim, ydim, zdim]):subv[i] for i in range(sbc)}

        endt = time.clock()

        print("done! ({0:.2f} kB in {1:.3f}s)\n".format(insize/1024.0, endt - startt))


    # find out how many molecules do we have
    # for each specie
    spnum = [0] * spc
    for i in subv:
        sp = subv[i]
        for j in range(spc):
            spnum[j] += sp[j]

    # print some info about the system
    print("      system shape       " \
          + str(xdim) + "x" + str(ydim) + "x" + str(zdim))
    print("      subvolumes count   " + str(sbc) + "\n")
    print("      species count      " + str(spc))
    print("      molecules count    " + str(sum(spnum)))
    for i in range(spc):
        print("        specie " + str(i) + "           " + str(spnum[i]))

    # write .vtk data to the output file
    print("\n    generating .vtk file.. ", end="");

    startt = time.clock()

    fpath = sys.argv[1][:-4]+".vtk" 
    f = open(fpath, "w+")

    writeheader(f, sum(spnum))
    for spi in range(spc):
        writepoints(f, subv, spi)
    f.write("\nPOINT_DATA " + str(sum(spnum)) + "\n")
    writecolors(f, spnum)
        
    f.close()
    outsize = os.path.getsize(fpath)

    endt = time.clock()
    print("done! ({0:.2f} kB in {1:.3f} s)".format(outsize/1024.0, endt - startt))
    
    
