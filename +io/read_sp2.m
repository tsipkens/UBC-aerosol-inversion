
% READ_SP2  Import SP2 data in the form of a list of timestamped particles.
%  
%  AUTHOR: Timothy Sipkens, 2024-07-16

function [mrbc, t] = read_sp2(fn)

% Read in SP2 data file.
disp('Reading SP2 data file...');

warning off;
in = readtable(fn);
warning on;

tools.textdone();

t = datetime(datestr(in.TimeStamp_sec_));  % corresponding time stamp

% mrbc = in.IncandMass_fg_;  % SP2 response (list of particles)
mrbc = 5e-7 .* in.IncandRelPeak .^ 0.9553;  % SP2 response (list of particles)

end
