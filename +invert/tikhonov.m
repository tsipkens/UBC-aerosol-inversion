
% TIKHONOV  Performs inversion using various order Tikhonov regularization. 
%  Regularization takes place in 2D. The type of regularization or prior
%  information added depends on the order of Tikhonov applied. 
% 
%  X = invert.tikhonov(A,B,LAMBDA) inverts the system A*X = B using 
%  Tikhonov regularization, where A is a model matrix or kernel and B is a 
%  data vector, and appling a regularization parameter LAMBDA. By default, 
%  the order of Tikhonov and thus the constructed Tikhonov matrix are given 
%  as the default of the invert.tikhonov_lpr(...) function. 
% 
%  X = invert.tikhonov(A,B,LAMBDA,L) uses the pre-computed Tikhonov matrix
%  L during the inversion. This matrix can be computed using the
%  invert.tikhonov_lpr(...) function. 
% 
%  X = invert.tikhonov(A,B,LAMBDA,ORDER,N) applies Tikhonov regularization
%  of the order specified by ORDER, where ORDER is an integer in the range 
%  [0,2]. This requires the quantity N, an integer representing the length 
%  of the first dimension of the solution, to build the Tikhonov matrix. 
%  This is relevant for full (not partial) grids. 
% 
%  X = invert.tikhonov(A,B,LAMBDA,ORDER,GRID) applies Tikhonov regularization
%  of the order specified by ORDER, as above. In the case of a partial
%  grid, an instance of the Grid class is required to build the Tikhonov
%  matrix in the place of an integer dimension.
% 
%  X = invert.tikhonov(...,BC) adds an input to control boundary
%  conditions.
% 
%  X = invert.tikhonov(...,XI) applies Tikhonov regularization using an
%  initial guess of XI. The default if not included is zeros. 
% 
%  X = invert.tikhonov(...,XI,SOLVER) applies Tikhonov regularization 
%  using the solver specified by the string in SOLVER. 
% 
%  [X,D] = invert.tikhonov(...) applies Tikhonov regularization, outputting
%  explicit inverse operator implied by this procedure, specifically such
%  that X = D*([B;0]). 
%  
%  [X,D,LPR0] = invert.tikhonov(...) applies Tikhonov regularization, also
%  outputting the unscaled Tikhonov matrix. 
% 
%  [X,D,LPR0,GPO_INV] = invert.tikhonov(...) applies Tikhonov
%  regularization, also outputting the inverse posterior covariance. 
% 
%  AUTHOR: Timothy Sipkens, 2018-11-21

function [x, D, Lpr0, Gpo_inv] = ... 
    tikhonov(A, b, lambda, order_L, n_grid, bc, xi, solver)  % BC before XI

x_length = size(A,2);

%-- Parse inputs ---------------------------------------------------------%
if ~exist('order_L','var'); order_L = []; end
    % if order not specified, use default of tikhonov_lpr

if ~exist('bc','var'); bc = []; end
if ~exist('xi','var'); xi = []; end % if initial guess is not specified
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%


%-- Get Tikhonov smoothing matrix ----------------------------------------%
if or(all(size(order_L)==[1,1]), iscell(order_L)) % if order is specified, build Lpr0
    Lpr0 = invert.tikhonov_lpr(...
        order_L, n_grid, x_length, bc);
else % is Lpr0 strucutre is provided directly
    Lpr0 = order_L;
end
Lpr = lambda .* Lpr0;


%-- Remove null data from data set ---------------------------------------%
%   This corresponds to data that does not impact reconstruction. 
keep = ~all(A == 0, 2);
A = A(keep, :);
b = b(keep);

%-- Choose and execute solver --------------------------------------------%
pr_length = size(Lpr0, 1);
[x,D] = invert.lsq(...
    [A;Lpr], [b;sparse(pr_length,1)], xi, solver);


%-- Uncertainty quantification -------------------------------------------%
if nargout>=4
    Gpo_inv = A'*A+Lpr'*Lpr;
end

end
%=========================================================================%

