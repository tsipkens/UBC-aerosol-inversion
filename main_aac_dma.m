
clear;
close all;
clc;

addpath tfer autils cmap;


%== (1) ==================================================================%
span = [20, 800; ...
    10, 600];  % span of grid
ne = [100, 125]; % number of elements in grid for each dimension

% Create an instance of Grid, with logarithmic spacing.
grid_x = Grid(span, ne, 'log');
grid_x.type = {'dm', 'da'};

ut_r = [2.3,2]; % point in line to cut upper triangle
ut_m = 1; % slope for line to cut upper triangle
lt_r = [1.5,2]; % point in line to cut lower triangle
lt_m = 1; % slope for line to cut upper triangle
grid_x = grid_x.partial(...
    ut_r, ut_m,...
    lt_r, lt_m); % convert to a partial grid


Gam = corr2cov(log10([1.5, 1.5]'), [1,0.98;0.98,1]);
phantom = Phantom('standard', [], log10([100, 120]), Gam);
x0 = phantom.eval(grid_x);  % evaluate phantom on grid_x

% Plot the phantom.
figure(1);
subplot(1, 2, 1);
grid_x.plot2d(x0);
colormap(viridis);  % apply primary colormap



%== (2) ==================================================================%
% Define a new grid for the measurements
span_b = span;
ne_b = [20, 65];
grid_b = Grid(span_b, ne_b, 'log');


% Get properites. 
prop_a = prop_aac();
prop_a = massmob.add(prop_a, 'zet', 3, 'rho100', 2160);

prop_a.Qa = prop_a.Qa(1) .* ones(size(grid_b.edges{1}));
prop_a.Qa(end-2:end) = 2.5e-5;

prop_a.Qexh = prop_a.Qexh(1) .* ones(size(grid_b.edges{1}));
prop_a.Qexh(end-2:end) = 2.5e-4;

prop_a.Qs = prop_a.Qa;
prop_a.Qsh = prop_a.Qexh;


prop_d = prop_dma();
prop_d.Qa = prop_d.Qa(1) .* ones(size(grid_b.edges{1}));
prop_d.Qa(end-2:end) = 2.5e-5;

prop_d.Qc = prop_d.Qc(1) .* ones(size(grid_b.edges{1}));
prop_d.Qc(end-2:end) = 2.5e-4;

prop_d.Qs = prop_d.Qa;
prop_d.Qm = prop_d.Qc;


% opts_a = struct();  opts_a.model = 'lt';

% Generate the kernel, use default CPMA properties. 
% A = kernel.build_grid(grid_b, grid_x, 1:3, ...
%     'charger', {}, 'dma', {prop_d}, 'aac', {prop_a});
A = kernel.build(grid_x, 1:4, 'charger', {}, ...
    'dma', {grid_b.elements(:,1), prop_d}, ...
    'aac', {grid_b.elements(:,2), prop_a});

A(A < 1e-3 .* max(max(A))) = 0;  % remove small values (speeds inversion)

figure(2);
grid_x.plot2d_marg(A(527,:)); % plot kernel for 527th data point
colormap(inferno);  % apply alternative colormap


%== (3) ==================================================================%
b0 = A * x0; % generate a set of data using the forward model

[b, Lb] = tools.get_noise(b0, 1e5); % corrupt data, assume peak counts ~1e5


%== (4) ==================================================================%
disp('Tikhonov inversion ...');
lambda = 0.2;  % regularization parameter
order = 1;  % order of Tikhonov matrix to be used
x_tk1 = invert.tikhonov(Lb*A, Lb*b, ...
    lambda, order, grid_x);  % Tikhonov solution
tools.textdone();
disp(' ');
disp(' ');
disp(' ');



%== (5) ==================================================================%
figure(1);
subplot(1, 2, 2);
grid_x.plot2d(x_tk1);  % plot Tikhonov solution
colormap(viridis);  % apply primary colormap

