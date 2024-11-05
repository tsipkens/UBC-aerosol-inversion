
% READ_CPOF  Read CPMA CPOF files. 
%  
%  AUTHOR: Timothy Sipkens, 2023-10-23
%  
%  TO ADD: Select a specific time interval in the file?

function [sp, prop, m_star, t_span] = read_cpof(fn, prop, off)

addpath tfer;

% Check whether to ignore when the SP2 is off. 
if ~exist('off', 'var'); off = []; end
if isempty(off); off = 1; end

if ~exist('prop', 'var'); prop = []; end
if isempty(prop)
    prop = prop_pma('cpma');
end

% Open CPMA CPOF file and read data.
opts = detectImportOptions(fn, ...
    'FileType', 'text');

warning off;
in = readtable(fn, opts);
warning on;

% Flag when the CPMA is off and remove data if relevant.
fl_off = in.Speed_rad_s_ < 100;
if off == 1
    in(fl_off,:) = [];  % remove off cases
end

% Use m_star to detect setpoint times. 
% u = unique([in.Mass_fg_, in.Rm], 'rows');
% m_star = u(:,1);
% Rm = u(:,2);
m_star = in.Mass_fg_(1);  Rm = in.Rm(1);  idx1 = 1;  idx2 = [];
for ii=2:height(in)
    if ~all([in.Mass_fg_(ii), in.Rm(ii)] == [m_star(end), Rm(end)])
        m_star = [m_star; in.Mass_fg_(ii)];
        Rm = [Rm; in.Rm(ii)];
        idx1 = [idx1; ii];
        idx2 = [idx2; ii-1];
    end
end
idx2 = [idx2; height(in)];
if idx1(2) > 15
    idx1(1) = 15;  % avoids some of initial ramp
else  % otherwise remove first "setpoint"
    idx1(1) = [];
    idx2(1) = [];
    m_star(1) = [];
    Rm(1)= [];
end


% Take flow rate as global median.
% This does not allow for flow rates to change in a given file.
Q = median(in.SampleFlow_lpm_);

% Update prop.
prop = prop_update_flow(prop, Q/1000/60);
prop.T = median(in.Temperature_C_) + 273.15; % average temperature
prop.p = median(in.Pressure_Pa_) ./ 101325; % average pressure



% Extract time span from CPOF file. 
% Assume only single time at each setpoint.
% Currently also do not account for buffer time between scans. 
t_span(length(m_star),2) = datetime();  % time span
idxs =(1:height(in))';
for ii=1:length(m_star)

    fl_ii = and(idxs >= idx1(ii), idxs <= idx2(ii));  % flag times for this sepoint
    
    % Get voltage over that interval (speed set below).
    % Then, remove outliers at beginning and end as identified for V.
    for jj=1:2
        V(ii) = median(in.Voltage_V_(fl_ii));
        [~,outl_V] = rmoutliers(in.Voltage_V_(fl_ii),'ThresholdFactor',1);
        fl_ii(fl_ii) = (~outl_V) .* fl_ii(fl_ii);
    end
    
    % Remove outliers at beginning and end as identified for omega.
    for jj=1:2
        omega(ii) = median(in.Speed_rad_s_(fl_ii));
        [~,outl_omega] = rmoutliers(in.Speed_rad_s_(fl_ii),'ThresholdFactor',0.5);
        fl_ii(outl_omega) = 0;
    end
    
    % Get time span after outlier removal. 
    idx1_post = find(fl_ii, 1, 'first');
    idx2_post = find(fl_ii, 1, 'last');

    t_span(ii,:) = [in.Date_Time(idx1_post),  in.Date_Time(idx2_post)];
end

% Add tfer_pma submodule. 
% Used to interpret setpoints below.
addpath('tfer/tfer-pma');
sp = get_setpoint(prop, 'V', V, 'omega', omega);


% Diagnostic plots.
% Uncomment to show plot.
%{
figure(gcf);
subplot(4, 1, 1);
plot(in.Date_Time, in.Voltage_V_);
yline(V);
set(gca, 'YScale', 'log');
ylim([1e1, inf]);
ylabel('Voltage');

subplot(4, 1, 2);
plot(in.Date_Time, in.Speed_rad_s_);
yline(omega);
set(gca, 'YScale', 'log');
ylim([1e2, inf]);
ylabel('Speed');

subplot(4, 1, 3);
plot(in.Date_Time, in.Mass_fg_);
hold on;
plot(t_span(:,1), [sp.m_star] .* 1e18, 'ro');
plot(t_span(:,2), [sp.m_star] .* 1e18, 'rx');
hold off;
yline(m_star);
xline(t_span(:,1), '--');
set(gca, 'YScale', 'log');
ylabel('m_star');

% Add shading for selected regions.
limy = ylim();
for ii=1:length(m_star)
    hold on;
    area([t_span(ii,1),t_span(ii,2)],[limy(2),limy(2)],'basevalue',limy(1), 'FaceColor', [0.9,0.9,0.9]);
    hold off;
end
h = get(gca,'Children');
set(gca,'Children',[h(end-length(m_star):end); h(1:end-length(m_star)-1)]);
ylim(limy);

subplot(4, 1, 4);
plot(in.Date_Time, in.Rm);
yline(Rm);
ylabel('Rm');
%}

end
