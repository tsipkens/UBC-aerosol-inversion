
% OVERLAY_PHANTOM  Overlay circles for std. dev. isolines for phantom (optionally, fits a unimodal phantom).
% Author: Timothy Sipkens, 2019-10-31
%=========================================================================%

function [pha, N, ci] = overlay_phantom(x_pha,grid,iso_levels,varargin)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('grid', 'var'); grid = []; end

if isa(x_pha, 'Phantom')
    pha = x_pha;
    if isempty(grid); grid = pha.grid; end
    N = [];
else
    [pha, N, ~, ci] = Phantom.fit(x_pha, grid);
end

if ~exist('iso_levels', 'var'); iso_levels = [1,2,3]; end  % plot 1, 2, and 3 sigma

if isempty(varargin); varargin = {'Color',[1,1,0,1]}; end
    % specify line properties (default: yellow, no transparency)
%-------------------------------------------------------------------------%


ylim_st = ylim;
xlim_st = xlim;

%-- Proceed with plotting ----------------%
mu = pha.mu; % local copies of parameters
Sigma = pha.Sigma;
p = pha.p;
for ll=1:pha.n_modes
    for jj=1:length(iso_levels)
        tools.overlay_ellipse(mu(ll,:),Sigma(:,:,ll),iso_levels(jj),varargin{:});
    end
    tools.overlay_line(grid,fliplr(mu(ll,:)),p(ll).Dm,varargin{:});
end


ylim(ylim_st); % restore original plot bounds
xlim(xlim_st);

if nargout==0; clear pha N; end % prevent unnecessary output

end

