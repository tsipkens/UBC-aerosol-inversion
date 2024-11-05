
% PLOT2D_SCATTER Scatter plot of data, where colour indicates counts.
% Author: Timothy Sipkens, 2019-11-28
%=========================================================================%

function [] = plot2d_scatter(vec1, vec2, b, cm)

% If not colormap specified, use grays.
if ~exist('cm', 'var'); cm = []; end
if isempty(cm); cm = gray(255); end


b = b ./ max(b); % scale data

% logb = log10(b); % compute log, allows for plotting of larger range
bmax = max(b);
bmin = min(b(~isinf(b)));
bmin = max(bmin,bmax-5); % at most, span four orders of magnitude
b = max(b,bmin);
bscl = max((b-bmin)/(bmax-bmin),bmin);

marker_size = 25 .* bscl + 0.2;
corder = 1 - bscl; % color is logscale

if ~exist('cm','var'); cm = []; end
if isempty(cm); cm = colormap('gray'); end

N = size(cm,1);
color = cm(round(corder.*(N-1)+1),:);
color(corder==1,:) = 1;

clf;
scatter(vec2, vec1,...
    marker_size, color, 'filled', 'MarkerEdgeColor', 'k');
    % marker_size is also logscale

set(gca,'XScale','log');
set(gca,'YScale','log');

%-- Add custom legend ----------------------%
hold on;
for ii=5:-1:0
    t0 = ii / 5.;
    h(6-ii) = scatter(NaN, NaN, 25 * t0 + 0.2, ...
        cm(round((1-t0) * (N-1) +1 ), :), 'filled', 'MarkerEdgeColor', 'k');
    text{6-ii} = num2str(t0);
end
% text{end} = ['<',text{end}];
hold off;

lgd = legend(h,text,'Location','northwest');
title(lgd,'{\times} b_{max}');
%-------------------------------------------%

end