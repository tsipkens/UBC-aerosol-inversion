
% Evaluates the JAS19 cases for a varying number of data points.

clear;
clc;
close all;


%-- Load colour schemes --------------------------------------------------%
addpath('cmap');
cm = load_cmap('YlGnBu',255);
cm_alt = cm;
load('inferno.mat');
cm = cm(40:end,:);
cm_b = cm;
load('viridis.mat');


%%
%-- Generate phantom (x_t) -----------------------------------------------%
%   High resolution version of the distribution to be projected to coarse
%   grid to generate x.
span_t = [10^-1.5,10^1.5;10,10^3]; % range of mobility and mass

phantom = Phantom('1',span_t);
x_t = phantom.x;
grid_t = phantom.grid;
nmax = max(x_t);
cmax = nmax;

figure(1);
phantom.plot;
colormap(gcf,[cm;1,1,1]);
caxis([0,cmax*(1+1/256)]);

hold on; % plots mg ridges of phantom
plot(log10(grid_t.edges{2}),...
    log10(phantom.mg_fun(grid_t.edges{2})),'w:');
hold off;

n_b_vec = [25,35,42,50,60,75];
for ii=1:length(n_b_vec)
    %%
    %-- Generate A matrix and b vector ---------------------------------------%
    n_b = [14,n_b_vec(ii)]; %[12,50]; %[17,35];
    span_b = grid_t.span;
    grid_b = Grid(span_b,...
        n_b,'logarithmic'); % should be uniform basis

    A_t = kernel.gen_grid(grid_b,grid_t,prop_pma,'Rm',3);
        % generate A matrix based on grid for x_t and b


    %%
    %--  Generate x vector on coarser grid -----------------------------------%
    n_x = [50,64]; % number of elements per dimension in x
        % [20,32]; % used for plotting projections of basis functions
        % [40,64]; % used in evaluating previous versions of regularization

    grid_x = Grid([grid_t.span],...
        n_x,'logarithmic');
    x0 = grid_x.project(grid_t,x_t); % project into basis for x

    disp('Transform to discretization in x...');
    B = grid_x.transform(grid_t); % evaluate matrix modifier to transform kernel
    A = A_t*B; % equivalent to integration, rebases kernel to grid for x (instead of x_t)
    A = sparse(A);
    disp('Complete.');
    disp(' ');

    figure(2);
    colormap(gcf,[cm;1,1,1]);
    grid_x.plot2d_marg(x0,grid_t,x_t);
    caxis([0,cmax*(1+1/256)]);


    %%
    %-- Generate data --------------------------------------------------------%
    b0 = A_t*x_t; % forward evaluate kernel

    %-- Corrupt data with noise ----------------------------------------------%
    b0(0<1e-10.*max(max(b0))) = 0; % zero very small values of b

    Ntot = 1e5;
    [b,Lb] = tools.get_noise(b0,Ntot);


    figure(5);
    colormap(gcf,cm_b);
    grid_b.plot2d_marg(b);

    figure(20);
    n2 = floor(grid_b.ne(1));
    n3 = floor(length(cm_b(:,1))/n2);
    cm_b_mod = cm_b(10:n3:end,:);
    set(gca,'ColorOrder',cm_b_mod,'NextPlot','replacechildren');
    b_plot_rs = reshape(b,grid_b.ne);
    semilogx(grid_b.edges{2},b_plot_rs.*Ntot);

    x_vec(ii).b = max(b);


    %-- Perform inversions -----------------------------------------------%
    run_inversions_a; % run regularization


    %-- Store values -----------------------------------------------------%
    x_vec(ii).tk0 = x_tk0;
    x_vec(ii).tk1 = x_tk1;
    x_vec(ii).tk2 = x_tk2;
    x_vec(ii).out_tk0 = out_tk0;
    x_vec(ii).out_tk1 = out_tk1;
    x_vec(ii).out_tk2 = out_tk2;
    x_vec(ii).out_mart = out_mart;
    x_vec(ii).out_two = out_two;
    x_vec(ii).out_two_mh = out_two_mh;
    x_vec(ii).chi_tk0 = chi.tk0;
    x_vec(ii).chi_tk1 = chi.tk1;
    x_vec(ii).chi_tk2 = chi.tk2;
    x_vec(ii).chi_mart = chi.mart;
    x_vec(ii).chi_two = chi.two;
    x_vec(ii).chi_two_mh = chi.two_mh;
    x_vec(ii).chi_lsq = chi.lsq;
end


%%
%-- Plot solution --------------------------------------------------------%
x_plot = x_vec(1).x_two;

figure(10);
colormap(gcf,[cm;1,1,1]);
grid_x.plot2d_marg(x_plot,grid_t,x_t);
caxis([0,cmax*(1+1/256)]);

%{
figure(13);
n1 = ceil(grid_x.ne(1)./20);
n2 = floor(grid_x.ne(1)/n1);
n3 = floor(240/n2);
cm_x = cm(10:n3:250,:);
set(gca,'ColorOrder',cm_x,'NextPlot','replacechildren');
x_plot_rs = reshape(x_plot,grid_x.ne);
semilogx(grid_x.edges{2},x_plot_rs(1:n1:end,:));

figure(10);
%}


%%
%-- Bar plot of results --------------------------------------------------%
figure(30);
chi_names = fieldnames(chi);
chi_vals = zeros(length(chi_names),1);
for ii=1:length(chi_names)
    chi_vals(ii) = chi.(chi_names{ii});
end

bar(chi_vals);
% ylim([0,20]);
% ylim([0,100]);
set(gca,'xticklabel',chi_names);


%%
%{
%-- Bar plot of times ----------------------------------------------------%
figure(40);
t_names = fieldnames(t);
t_vals = zeros(length(t_names),1);
for ii=1:length(t_names)
    t_vals(ii) = mean(t.(t_names{ii}),2);
end

bar(t_vals);
set(gca,'xticklabel',t_names);
set(gca,'yscale','log');



%%
%-- Plot marginal distributions ------------------------------------------%
figure(31);
clf;
dim = 2;

grid_t.plot_marginal(x_t,dim);
grid_x.plot_marginal(...
    {x_Tk1,x_init,x_MART,x_Two,x_TwoMH},dim,x0);



%%
%-- Plot conditional distributions ---------------------------------------%
figure(31);
clf;
dim = 2;
ind_plot = 25;

grid_x.plot_conditional(...
    {x0,x_Tk1,x_init,x_MART,x_Two,x_TwoMH},dim,ind_plot,x0);
%}
