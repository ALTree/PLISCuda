## usage: [data, times] = parseserie(statedir, specie, sbc)
##
## 'statedir' is the path of the dir with the .dat filesep
## 'specie' is the specie of interest (i.e. its index)
## 'sbc' is the number of subvolumes of the model
##
## parseserie returns a data matrix and a vector of times.
##
## The n-th row of the matrix contains the molcount (vs time) of the specie 
## 'specie', in the n-th subvolume.
##
## 'times' is a vector with the times of the snapshots.
##
## plot(times, data(42,:)) can be used to plot the number of mols of specie
## 'specie' in subvolume 42 over time.

function [data, times] = parseserie (dirpath, specie, sbc)

datfiles = glob(strcat(dirpath, "/sim*"));
times = zeros(1, length(datfiles));

# We need a first pass on the files because the timestamp is in the
# filename and we don't want to assume that the OS will give us the files
# in timestamp-order (in fact, it won't).
#
# Parse every file and put the timestamps in the 'times' vector. 
for i = 1:length(datfiles)
  [_, t] = parsestate(datfiles{i});
  times(i) = t;
end

# and now we know how long the timeserie is, and what the timestep is.

step = floor(min(times));  # this is also the smallest time
maxt = floor(max(times));  # biggest time


# preallocate data matrix and fill it

data = zeros(sbc, floor(maxt/step));
for i = 1:length(datfiles)
  [state, t] = parsestate(datfiles{i});
  data(:, floor(t/step)) = state(:,specie);
end

# sort times or plot(times, ..) won't work
times = floor(sort(times));

end