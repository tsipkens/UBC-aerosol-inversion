
% X2Y  A generic function to transform a distribution to a different space.
%  Examples include converting mass-mobility distributions to
%  effective-density mobility distributions or converting mrBC-mp
%  distributions to frBC-mp distributions. 
%  
%  Y = tools.x2y(X, GRID_X) converts a mass-mobility distribtuion to
%  an effective density-mobility distribution. By default DIM = 2, SPAN_Y
%  is estimated from the input data and grid span, and N_Y = 600. See notes
%  below for variable definitions.
%  
%  Y = tools.x2y(X, GRID_X, FUN) add an input for a transformation 
%  function_handle. Expected format is `FUN = @(A,B)...` where A is the
%  quantity in the first dimension of the grid and B is the quantity in the
%  second dimension. FUN must be a linear function in logspace for A and B. 
%  
%    FUN examples:
%    PMA-DMA:
%       A. mass-mobility > effective density-mobility
%           FUN = @(a, b) 6 .* a ./ (pi .* b .^ 3) .* 1e9;
% 
%    PMA-SP2
%       B. mrBC-mp > frBC-mp
%          FUN = @(a, b) a ./ b;
%  
%  Y = tools.x2y(X, GRID_X, FUN, DIM) applies the transformation while
%  preserving the quantity in the DIM dimension of the input grid. For
%  example, for mass-mobility distributions and the default of DIM = 2, the
%  transformation will be applied while preserving the mobility diameter,
%  which will remain as the second dimension in the output grid.
%  
%  Y = tools.x2y(..., SPAN_Y) explicitly states the span of the transformed
%  quantity on the output grid. By default, the span is estimated based on
%  the max and min of x values that are within three order-of-magnitude of
%  max{x}.
%  
%  Y = tools.x2y(..., N_Y) explicitly states the number of elements
%  in the effective density dimension. 
%  
%  [Y, GRID_Y] = tools.x2y(...) outputs the grid on which the
%  effective density-mobility distribution is defined.
%  
%  [Y, GRID_Y, T] = tools.x2y(...) add output for transformation matrix.
%  
%  ------------------------------------------------------------------------
%  
%  NOTE: This function requires that FUN be linear in logspace (e.g.,
%  allowing for power laws) and, by extension, to be monotonic.
%  
%  AUTHOR: Timothy Sipkens, 2019-05-17

function [y, grid_y, T] = x2y(x, grid_x, fun, dim, span_y, n_y)

%-- Parse inputs -----------------------------------%
if ~exist('fun', 'var'); fun = []; end
if isempty(fun)  % by default use mass-mobility distribution
    fun = @(a, b) 6 .* a ./ (pi .* b .^ 3) .* 1e9;
end

% By default, set first dimension as common dimension.
if ~exist('dim', 'var'); dim = []; end
if isempty(dim); dim = 2; end
dim2 = 3 - dim;  % other dimension, i.e., 2 -> 1 or 1 -> 2

if ~exist('n_y', 'var'); n_y = []; end
if isempty(n_y); n_y = 600; end  % can be large as conversion is simple

if ~exist('span_y', 'var'); span_y = []; end
if isempty(span_y)  % estimate from existing spans and x values
    f0 = fun(grid_x.elements(:, 1), grid_x.elements(:, 2));  % compute fun for all elements
    
    % f_sig = x > max(x) ./ 1e4;  % flag significant x values

    f2 = log10(max(f0(f_sig)));  % upper bound
    f2 = ceil(10 .^ (f2 - floor(f2)) .* 10) ./ 10 .* ...  % pre-factor
        10 .^ floor(f2);  % exponent
    
    f1 = log10(min(f0(f_sig)));  % lower bound
    if f1 == -Inf; f1 = log10(f2 ./ 1e3); end  % if zero, span 5 orders of magnitude instead
    f1 = floor(10 .^ (f1 - floor(f1)) .* 10) ./ 10 .* ...  % pre-factor
        10 .^ floor(f1);  % exponent
    
    span_y = [f1, f2];
end
%---------------------------------------------------%


% Quantities relevant for new grid.
y_min = span_y(1);  % get span for neww quantity
y_max = span_y(2);
y_n = logspace(log10(y_min), ...
               log10(y_max), n_y);  % discretized y space

% Generate new grid for transformed space.
% Generate one new dimension while keeping dimension DIM.
if dim == 2
    grid_y = Grid([y_min,y_max; grid_x.span(2,:)],...
        [n_y, length(grid_x.edges{2})], 'logarithmic');
else
    grid_y = Grid([grid_x.span(1,:); y_min,y_max],...
        [length(grid_x.edges{1}), n_y], 'logarithmic');
end

x_rs = grid_x.reshape(x);

% Transpose such that second dimension is always the one preserved.
% Simplifies some of the below procedure.
if dim == 1; x_rs = x_rs'; end


%== Loop over mobility diameter ==========================================%
%   (i.e., consider conditional mass distributions)
n_dim = grid_x.ne(dim);  % number of elements for DIM of grid_x
y = zeros(grid_x.ne(dim), length(y_n));  % intialize transforms distribution

for ii = 1:n_dim
    
    % Initialize transformation kernel of appropriate size.
    T = zeros(grid_y.ne(dim2), grid_x.ne(dim2));
    
    % Convert x nodes to effective density for iith mobility.
    if dim == 2
        y_old = fun(grid_x.nodes{1}, grid_x.edges{2}(ii));
    else
        y_old = fun(grid_x.edges{1}(ii), grid_x.nodes{2});
    end
    
    % Reverse order if y_old is decreasing.
    % Overlapping elements considerations below requires 
    % monotonic increase.
    f_reverse = 0;
    if y_old(2) < y_old(1)
        y_old = fliplr(y_old);
        f_reverse = 1;
    end
    
    % Take logarithm.
    y_old = max(y_old, 0);  % remove negative entries (not allowed on log grid)
    y_old = log10(y_old);
    
    % Loop over masses computing overlap between new and old elements.
    for jj = 1:grid_x.ne(dim2)
        T(:,jj) = max(...
            min(log10(grid_y.nodes{dim2}(2:end)), y_old(jj+1)) - ... % upper bound
            max(log10(grid_y.nodes{dim2}(1:(end-1))), y_old(jj)) ... % lower bound
            , 0) ./ ...
            (log10(grid_y.nodes{dim2}(2:end)) - ...
            log10(grid_y.nodes{dim2}(1:(end-1)))); % normalize by y bin size
    end
    
    % Re-flip T if f_reverse flagged above.
    if f_reverse; T = fliplr(T); end
    
    % Multiply transformation by mass-mobility distr.
    y(ii, :) = T * x_rs(:, ii);
end
%=========================================================================%


% Format data for output
if dim == 2; y = y'; end  % transpose for DIM = 2 to preserve original dimension
y = y(:);

end

