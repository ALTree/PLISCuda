## plotstate2.m

## Author:  <a.donizetti@cineca.it>
## Created: 2016-07-26

function [] = plotstate2 (state, dims)

xdim = dims(1);
ydim = dims(2);
zdim = dims(3);

binsz = zeros(ydim, xdim); % flatten along z dimension
binsy = zeros(zdim, xdim); % flatten along y dimension
binsx = zeros(ydim, zdim); % flatten along x dimension

for i=1:length(state)
	[x, y, z] = ind2sub([xdim, ydim, zdim], i);
	binsz(y, x) = binsz(y, x) + state(i); 
	binsy(z, x) = binsy(z, x) + state(i);
	binsx(y, z) = binsx(y, z) + state(i);
end

subplot(3, 2, [1 2]);
imagesc(binsz);
colorbar;
title("flatten z");

subplot(3, 2, [3 4]);
imagesc(binsy);
colorbar;
title("flatten y");

subplot(3, 2, 5);
imagesc(binsx, [0 max(max(binsx))]);
axis([0 ydim 0 zdim], "equal");
colorbar;
title("flatten x");

endfunction

