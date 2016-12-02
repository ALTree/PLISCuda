## utils.m

## countmols

## usage: mols = countmols(state, dims, low, high)
##
## Returns the number of mols in 'state' that are included in the box delimited
## by subvolumes low = [xl, yl, zl] and high = [xh, yh, zh] (extremes included).
##
## Example (for a system with a 10x10x10 topology):
##   mols = countmols(state, [10 10 10], [0 0 0], [4 4 4])
## Will return an array with the number of mols of each specie that are included
## in the box delimited by the the point (0,0,0) and the point (5,5,5).

function [mols] = countmols (state, dims, low, high)

[xdim, ydim, zdim] = num2cell(dims){:};
[xl, yl, zl] = num2cell(low){:};
[xh, yh, zh] = num2cell(high){:};

spc = length(state(1,:));  # species count
sbc = length(state(:,1));  # subvolumes count

mols = zeros(1, spc);

for i = 1:sbc
  [x, y, z] = ind2sub(dims, i);
  x = x-1; y = y-1; z = z-1;
  if (
    (x >= xl && x <= xh) &&
    (y >= yl && y <= yh) &&
    (z >= zl && z <= zh))
    mols = mols .+ state(i,:);
    end
end



end