## Copyright (C) 2016 
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


## usage: [state, t] = parsestate(filename)
##
## Given a system state file with N species and V subvolumes, returns an NxV
## matrix with species on the columns and the subvolumes on the rows (in state)
## and the corresponding simulation time (in t).


## Author:  <a.donizetti@cineca.it>
## Created: 2016-07-26

function [state, t] = parsestate (filename)

fid = fopen(filename, "r");

sbc = str2num(strsplit(fgets(fid), "="){2});
spc = str2num(strsplit(fgets(fid), "="){2});
t = str2double(strsplit(fgets(fid), "="){2});
fgetl(fid); % eat newline

fmtstr = "";
for spi = 1:spc
  fmtstr = strcat(fmtstr, " %d");  # this is incredibly inefficient
end

state = zeros(sbc, spc);

for i = 0:(sbc-1)
  v = fscanf(fid, fmtstr, spc)'; 
  state(i+1, :) = v';
end

fclose(fid);

endfunction















