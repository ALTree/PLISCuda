import sys


def fill_line(neigh, inds):
    # fill first and last element
    first, last = inds[0], inds[-1]
    neigh[first] = neigh[first].union({inds[1]})
    neigh[last] = neigh[last].union({inds[-2]})

    # fill the middle ones
    for i in range(1, len(inds)-1):
        neigh[inds[i]] = neigh[inds[i]].union({inds[i-1], inds[i+1]})

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("usage: $ python3 genTopology.py fileName xdim [ydim]")
        sys.exit(0)
            
    dims = [int(n) for n in sys.argv[2:]]
    if len(dims) == 1: # if we only got one dimension, fill with ydim = 1
        dims.append(1)

    xdim, ydim = dims[0:2]
    # initialize neigh dict so we don't have to put
    #     if i in neigh:
    # everywhere when we update it
    neigh = {i:set() for i in range(0, xdim*ydim)}

    # fill neigh by lines
    for y in range(0, ydim):
        fill_line(neigh, range(xdim*y, xdim*(y+1)))

    # fill neigh by colums
    if ydim > 1:
        for x in range(0, xdim):
            fill_line(neigh, range(x, x + (ydim*xdim), ydim))

    # sanity checks
    lens = [len(neigh[i]) for i in neigh]

    if dims[1] == 1:
        assert lens.count(1) == 2
        assert lens.count(2) == xdim-2
        assert lens.count(3) == 0
        assert lens.count(4) == 0
    else:
        assert lens.count(1) == 0
        assert lens.count(2) == 4
        assert lens.count(3) == 2*(xdim + ydim - 4)
        assert lens.count(4) == (xdim*ydim) - (2*(xdim + ydim - 4) + 4)

    # write the neigh dict to file
    with open(sys.argv[1], 'w') as f:
        for k in neigh:
            f.write(str(k) + ": ")
            for v in neigh[k]:
                f.write(str(v) + " ")
            f.write("\n")


