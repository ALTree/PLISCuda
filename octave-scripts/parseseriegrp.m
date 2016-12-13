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

function [data, times] = parseseriegrp (dirpath, specie, sbc, dims, gp)

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

data = zeros(dims(1)/gp, floor(maxt/step));
counter = 0;

for i = 1:length(datfiles)
  counter = counter + 1;
  
  [state, t] = parsestate(datfiles{i});
  
  for j = 0:gp:(dims(1)-1)
    mols = countmols(state, dims, [j 0 0], [j + gp - 1, dims(2)-1, dims(3)-1]);
    data((j/gp)+1, floor(t/step)) = mols(specie);
  end
  
  printf("Processed %d/%d\n", counter, length(datfiles));
  fflush(stdout);
  
end

# sort times or plot(times, ..) won't work
times = floor(sort(times));

end