## Copyright (C) 2016 Alberto Donizetti
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## usage: genstate(dims, mols, low, high)
##
## Generates a state file for a model with a dims = [x, y, z] topology.
## More specifically, genstate, for every subvolume between 
##   low = [xl, yl, zl] 
## and 
##   high = [xh, yh, zh]
## (extremes included) sets the state to mols = [specie1 specie2 ...].
## Every other subvolume gets a [0 0 ...] state.

## Author: Alberto Donizetti <adonizetti@nadonizet230119>
## Created: 2016-11-29

function [] = genstate (dims, mols, low, high)

# extract singular dimension from input arrays

[xdim, ydim, zdim] = num2cell(dims){:};
[xl, yl, zl] = num2cell(low){:};
[xh, yh, zh] = num2cell(high){:};

subvc = xdim*ydim*zdim;


# prepare fprintf format strings

fmtstr = "%d:";  # fmtstr for subv where we'll write 'mols'
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


# generate file header

fid = fopen("state.txt", "w");
fprintf(fid, "subvolumes: %d\n", subvc);
fprintf(fid, "species: %d\n\n", length(mols));


# generate subvolumes state lines

for i = 0:subvc-1
  [x, y, z] = ind2sub(dims, i+1);
  x = x-1; y = y-1; z = z-1;
  if (
    (x >= xl && x <= xh) &&
    (y >= yl && y <= yh) &&
    (z >= zl && z <= zh))
    fprintf(fid, fmtstr, i, mols);
    else
      fprintf(fid, zstr, i, 2);
    end
end


fclose(fid);

endfunction
