## plostate3.m

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