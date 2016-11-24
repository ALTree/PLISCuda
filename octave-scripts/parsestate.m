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


## Author:  <a.donizetti@cineca.it>
## Created: 2016-07-26

function [state] = parsestate (filename)

fid = fopen(filename, "r");

sbc = str2num(strsplit(fgets(fid), "="){2});
spc = str2num(strsplit(fgets(fid), "="){2});
time = str2double(strsplit(fgets(fid), "="){2});
fgetl(fid); % eat newline

fmtstr = "";
for spi = 1:spc
  fmtstr = strcat(fmtstr, " %d");
end

state = zeros(sbc, spc);

for i = 0:(sbc-1)
  v = fscanf(fid, fmtstr, spc)'; 
  state(i+1, :) = v';
end

endfunction















