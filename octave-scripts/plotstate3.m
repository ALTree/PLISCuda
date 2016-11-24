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

function [] = plotstate3(state, dims)

xdim = dims(1);
ydim = dims(2);
zdim = dims(3);

xv = zeros(1, 100);
yv = zeros(1, 100);
zv = zeros(1, 100);

for i=1:length(state)
	[x, y, z] = ind2sub([xdim, ydim, zdim], i);
	if state(i) != 0
		xv = [xv, (x-1)];
		yv = [yv, (y-1)];
		zv = [zv, (z-1)];
	end
end

scatter3(xv, yv, zv, 1, 'b');
axis([0 205 0 5 0 5], "equal");

endfunction