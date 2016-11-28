folder = "~/Downloads/out6000/";
chdir(folder);

specie = 9;  # specie of interest
subvc = 60;  # number of subvolumes

datfiles = glob("sim*")
times = zeros(1, length(datfiles));

for i = 1:length(datfiles)
  [_, t] = parsestate(datfiles{i});
  times(i) = t;
end

step = floor(min(times));  # this is also the smallest time
maxt = floor(max(times));  # biggest time

data = zeros(subvc, maxt/step);

for i = 1:length(datfiles)
  [state, t] = parsestate(datfiles{i});
  data(:, floor(t/step)) = state(:,specie);
end

times = floor(sort(times));

clear i maxt step specie subvc t datfiles _ folder state