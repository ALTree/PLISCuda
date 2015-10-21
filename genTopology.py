import sys

def add_neighs(sub):
    """
    add_neighs addd neighbours to the subvolume sub. 
    Ignores borders issues, we'll filter out invalid 
    neighbours later.
    """
    global topology

    for i in [-1, +1]:
        x, y, z = sub
        topology[sub].add((x+i, y, z))
        topology[sub].add((x, y+i, z))
        topology[sub].add((x, y, z+i))

def is_valid(neigh):
    """
    is_valid returns True iff neigh is a valid subvolume
    for the topology we are working in.
    """
    global xdim, ydim, zdim

    if -1 in neigh:
        return False
    
    x, y, z = neigh
    return (x < xdim and y < ydim and z < zdim)

def linearize(sub):
    """
    linearize computes a linear ID from an (x,y,z) subvolume ID
    """
    global xdim, ydim, zdim

    x, y, z = sub
    return x + y*xdim + z*xdim*ydim
        
if len(sys.argv) < 5:
    print("usage: python3 $genTopology.py x y z filename.txt")
    sys.exit()

xdim, ydim, zdim = [int(n) for n in sys.argv[1:-1]]

topology = {(x, y, z):set()
            for x in range(0,xdim)
            for y in range(0,ydim)
            for z in range(0,zdim)}

for i in topology:
    add_neighs(i)

# filter invalid neighbours
for i in topology:
    topology[i] = set(filter(is_valid, topology[i]))

# linearize keys
for i in topology.copy():
    topology[linearize(i)] = topology[i]
    topology.pop(i)

# linearize values
for i in topology.copy():
    topology[i] = {linearize(sub) for sub in topology[i]}

# sanity checks
lens = [len(x) for x in topology.values()]

if xdim > 1 and ydim == 1 and xdim == 1:
    assert lens.count(1) == 2
    assert lens.count(2) == xdim-2
    assert lens.count(3) == 0
    assert lens.count(4) == 0

if xdim > 1 and ydim > 1 and zdim == 1:
    assert lens.count(1) == 0
    assert lens.count(2) == 4
    assert lens.count(3) == 2*(xdim + ydim - 4)
    assert lens.count(4) == (xdim*ydim) - (2*(xdim + ydim - 4) + 4)

with open(sys.argv[4], 'w') as f:
    for i in sorted(topology):
        f.write(str(i) + ": ")
        for v in topology[i]:
            f.write(str(v) + " ")
        f.write("\n")

