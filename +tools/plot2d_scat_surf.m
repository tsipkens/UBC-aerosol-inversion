
% PLOT2D_SCAT_SURF  Takes a scattering of points and estimates a surfaces. 
%  Adds interpolation to fill in gaps, which limits accuracy. 
%  
%  Expects log-axes.
%  
%  AUTHOR: Timothy Sipkens, 2024-05-14

function plot2d_scat_surf(x, y, z)

F = scatteredInterpolant(x, y, z, 'natural', 'nearest');

x1 = logspace(log10(min(x)), log10(max(x)), 100);
y1 = logspace(log10(min(y)), log10(max(y)), 101);
[xq, yq] = meshgrid(x1, y1);

zq = F(xq, yq);

imagesc(x1, y1, zq);
hold on;
s = scatter(x, y, 'filled', 'w', 'SizeData', 4);
alpha(s, .5);
hold off;
set(gca, 'XScale', 'log', 'YScale', 'log');

xlim([min(x), max(x)]);
ylim([min(y), max(y)]);

end