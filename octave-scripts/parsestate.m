## usage: [state, t] = parsestate(filename)
##
## Given a system state file with N species and V subvolumes, returns an NxV
## matrix with species on the columns and subvolumes on the rows (in state)
## (and the corresponding simulation time, in t).

function [state, t] = parsestate (filename)

fid = fopen(filename, "r");

sbc = str2num(strsplit(fgets(fid), "="){2});
spc = str2num(strsplit(fgets(fid), "="){2});
t = str2double(strsplit(fgets(fid), "="){2});
fgetl(fid); # eat newline

fmtstr = "";
for spi = 1:spc
  fmtstr = strcat(fmtstr, " %d");  # this is incredibly inefficient
end

state = zeros(sbc, spc);

for i = 1:sbc
  state(i, :) = fscanf(fid, fmtstr, spc);
end

fclose(fid);

endfunction















