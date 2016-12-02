## genstate.m

## usage: genstate(dims, mols, low, high)
##
## Generates a state file for a model with a dims = [x, y, z] topology.
## More specifically, genstate, for every subvolume between 
##   low = [xl, yl, zl] 
## and 
##   high = [xh, yh, zh]
## (extremes included) sets the state to mols = [specie1 specie2 ...].
## Every other subvolume gets a [0 0 ...] state.
## Indices are 0-based.

## Author: Alberto Donizetti <adonizetti@nadonizet230119>
## Created: 2016-11-29

function [] = genstate (dims, mols, low, high)

[xdim, ydim, zdim] = num2cell(dims){:};
[xl, yl, zl] = num2cell(low){:};
[xh, yh, zh] = num2cell(high){:};

subvc = xdim*ydim*zdim;


# prepare fprintf format strings

fmtstr = "%d:";  # fmtstr for subv where we'll actually write 'mols'
for i = 1:length(mols)
  fmtstr = strcat(fmtstr, " %d");
end
fmtstr = strcat(fmtstr, "\n");

zstr = "%d:";  # fmtstr for subv where we'll only write zeros
for i = 1:length(mols)
  if i == 12
    zstr = strcat(zstr, " %d");
  else
    zstr = strcat(zstr, " 0");
  end
end
zstr = strcat(zstr, "\n");


# write file header

fid = fopen("state.txt", "w");
fprintf(fid, "subvolumes: %d\n", subvc);
fprintf(fid, "species: %d\n\n", length(mols));


# write subvolumes state lines

for i = 0:subvc-1
  [x, y, z] = ind2sub(dims, i+1);
  x = x-1; y = y-1; z = z-1;
  if (
    (x >= xl && x <= xh) &&
    (y >= yl && y <= yh) &&
    (z >= zl && z <= zh))
    fprintf(fid, fmtstr, i, mols);
    else
      fprintf(fid, zstr, i, 25);
    end
end


fclose(fid);

endfunction
