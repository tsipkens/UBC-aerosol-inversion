
% READ_CPMA  Read a CPMA file. 
%  
%  AUTHOR: Timothy Sipkens, 2024-07-15

function [sp, prop, t, ave] = read_cpma(fn, prop)

disp('Reading CPMA file ...');

addpath('tfer\tfer-pma');

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)
    prop = prop_pma('cpma');
end

% Open CPMA CPOF file and read data.
opts = detectImportOptions(fn, ...
    'FileType', 'text', 'VariableNamesLine', 8, ...
    'Delimiter', {'\t'}, 'Whitespace', '\b');
opts.DataLines = [9, Inf];

warning off;
in = readtable(fn, opts);
warning on;

% Make sure Datum_ is a number (can be erroneous from readtable).
if strcmp(opts.VariableOptions(1).Type, 'char')
    in.Datum_ = str2double(in.Datum_);
end
in(isnan(in.Datum_), :) = [];

%-- Read header lines -----------------%
fid = fopen(fn, 'r');
head = textscan(fid, '%s', 'headerlines', 0);
fclose(fid);

idx = find(strcmp('Density', head{1}));
rho0 = str2double(head{1}{idx(1) + 2}) ./ 1e3;

idx = find(strcmp('Fractal', head{1}));
zet = str2double(head{1}{idx(1) + 2});

prop = massmob.add(prop, 'zet', zet, 'm100', rho0 .* 1e-18);

idx = find(strcmp('(lpm):', head{1}));
Q = str2double(head{1}{idx(1) + 1});
prop = prop_update_flow(prop, Q/1000/60);

idx = find(strcmp('Temperature', head{1}));
prop.T = str2double(head{1}{idx(1) + 2});

idx = find(strcmp('Pressure', head{1}));
prop.p = str2double(head{1}{idx(1) + 2}) ./ 101325;

idx = find(strcmp('Averaging', head{1}));
ave = str2double(head{1}{idx(1) + 3});
%--------------------------------------%


sp = get_setpoint(prop, 'omega', in.Speed_rad_s_, 'V', in.Voltage_V_);

t = in.Time;

tools.textdone();

end
