
clear;
clc;
close all;


%-- Load colour maps -----------------------------------------------------%
addpath cmap tfer tfer/tfer-pma;
cm_alt = bupu(255);
cm_b = inferno(255);
cm_b = cm_b(40:end,:);
cm_div = rdbu(200);
cm = magma;
cm = flipud(cm);
cm = [1,1,1;cm];



%%
%== (1) ==================================================================%
%   Generate phantom (x_t).
%   High resolution version of the distribution to be projected to coarse
%   grid to generate x.
span_t = [ ...
    10^-2,10^2; ...  % range of mobilities
    10^-2,10^2]; % range of masses

grid_t = Grid(span_t,...
    [540,540],'logarithmic');
grid_t = grid_t.partial(0,1);
phantom = Phantom('distr-sp2-2', grid_t);
x_t = phantom.x;
nmax = max(x_t);
cmax = nmax;

%--  Generate x vector on coarser grid -----------------------------------%
n_x = [72,72]; % number of elements per dimension in x

grid_x = Grid([grid_t.span],...
    n_x,'logarithmic');
grid_x = grid_x.partial(0,1);
x0 = grid_x.project(grid_t,x_t); % project into basis for x

% phantom_a = Phantom('distr-sp2',grid_x);
% x_a = phantom_a.x;

figure(1);
tools.plot2d_marg_b(phantom.grid,phantom.x)
colormap(gcf,cm);
caxis([0,cmax*(1+1/256)]);

hold on; % plots mg ridges of phantom
plot(log10(grid_t.edges{2}),...
    log10(phantom.mg_fun(grid_t.edges{2})),'w:');
hold off;



%%
%== (2) ==================================================================%
%   Generate A matrix. Note that here a dense kernel is computed 
%   for data synthesis in Step 3. 

n_b = [60,15]; %[14,50]; %[17,35];
span_b = grid_t.span;
grid_b = Grid(span_b,...
    n_b,'logarithmic'); % grid for data
% grid_b = grid_b.partial(0,1);

prop = prop_pma;
[A_t,sp] = kernel.gen_pma_sp2_grid(grid_b, grid_t, prop, 'Rm', 3);
    % generate A matrix based on grid for x_t and b

disp('Transform to discretization in x ...');
B = grid_x.transform(grid_t); % evaluate matrix modifier to transform kernel
A = A_t*B; % equivalent to integration, rebases kernel to grid for x (instead of x_t)
A = sparse(A);
tools.textdone();
disp(' ');

figure(2);
colormap(gcf,cm);
d0 = grid_x.marginalize(x0);
tools.plot2d_marg_b(grid_x,x0,d0{1},grid_t,x_t);
caxis([0,cmax*(1+1/256)]);



%%
%== (3) ==================================================================%
%   Generate data.

b0 = A_t * x_t; % forward evaluate kernel

%-- Corrupt data with noise ----------------------------------------------%
b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

% Note: The total number of particles, Ntot, scales the pdf to equal the 
% number distribution, d2N/dlogA*dlogB. It is also equal to the 
% concentration of particles reduced by the product of the flow 
% rate and overall collection time, that is Ntot = N*Q*t. 
Ntot = 1e5; % converts pdf to counts
[b,Lb] = tools.get_noise(b0,Ntot);

figure(5);
tools.plot2d_scatter(...
    grid_b.elements(:,1),grid_b.elements(:,2),b,cm_b);

figure(20);
grid_b.plot2d_sweep(Ntot.*b,cm_b,2);
xlabel('{{\itm}_{p} [fg]}');
ylabel('{d{\itN}/dlog {\itm}_{p}}');



%%
%== (4) ==================================================================%
%   Invert.

s1 = 0.25;
s2 = s1;
R12 = 0.99;
Gd = [s1^2,R12*s1*s2; R12*s1*s2,s2^2];
lambda_ed = 1.2e0;

disp('Exponential distance ...');
Lpr0 = invert.exp_dist_lpr(Gd,grid_x.elements(:,2),grid_x.elements(:,1));
[x_ed,~,Lpr0] = invert.exp_dist(...
    Lb*A,Lb*b,lambda_ed,Gd,...
    grid_x.elements(:,2),grid_x.elements(:,1));
% [x_ed,lambda_ed,out] = optimize.exp_dist_op(...
%         Lb*A,Lb*b,[1e1,1e4],Gd,...
%         grid_x.elements(:,2),grid_x.elements(:,1));
tools.textdone();
disp(' ');



%%
%== (5) ==================================================================%
%   Visualize the results.

x_plot = x_ed;

%-- Plot retrieved solution --------------%
figure(10);
colormap(gcf,cm);
tools.plot2d_marg_b(grid_x,Ntot.*x_plot,[],grid_t,x_t);
tools.overlay_line(grid_x,[0,0],1,'k--');
% colorbar;
caxis([0,Ntot*cmax*(1+1/256)]);
colormap(cm);


