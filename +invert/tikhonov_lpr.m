
% TIKHONOV_LPR  Generates Tikhonov smoothing operators/matrix, L. 
% 
%  NOTE: This is a companions function to invert.tikhonov and does not 
%  performinversion. Further, this function does not include the effect of  
%  the regularization parameter.
%  
%  LPR0 = invert.tikhonov_lpr(ORDER, N, X_LENGTH) computes the  
%  bidimensional Tikhonov matrix with smoothing of order, ORDER, for a 2D
%  grid having N elements in the first dimension and a X_LENGTH total
%  elements. Values for ORDER at there simplest are integer values
%  corresponding to the type of derivative to be minimized. For example,
%  ORDER = 1 creates a matrix that estimates the first derivative of some
%  X. 
%  
%  LPR0 = invert.tikhonov_lpr(ORDER, GRID, []) replaces N and X_LENGTH with
%  an instance of the Grid or PartialGrid class, extracting the necessary
%  information and, in some instances, allowing for Tikhonov matrix for
%  grids that have missing elements (i.e., a PartialGrid). 
%  
%  LPR0 = invert.tikhonov_lpr(..., BC) adds a variable to specify the type
%  of boundary condition for certain methods. For example, BC = 0 pushed
%  the Tikhonov reconstruction to be zero at the boundaries. The exact
%  effect of BC depends on the ORDER and the specific evaluation method. In
%  general, the default is to match the boundary condition to the order,
%  i.e., BC = ORDER. 
%  
%  [LPR0, LPR1, LPR2] = invert.tikhonov(...) adds outputs corresponding to
%  subcomponent of the Tikhonov matrix (e.g., the differential operators in
%  different dimensions). 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2020-02-05
%  UPDATE: Timothy Sipkens, 2024-03-31

function [Lpr0, Lpr1, Lpr2] = tikhonov_lpr(order, n_grid, x_length, bc)

if ~exist('order', 'var'); order = []; end
if isempty(order); order = 1; end

% Parse extra variants, if relevant.
if iscell(order)
    var = order(2:end);
    order = order{1};
else; var = {0};
end

if ~exist('bc', 'var'); bc = []; end
if isempty(bc); bc = floor(order); end  % by default, match order

if ~isa(n_grid, 'Grid')
    n = n_grid;
    if and(mod(x_length, n_grid)~=0, order~=0) % error if dimensions don't make sense
        error('Error: x_length must be integer multiple of n.');
    end
end


%-- Generate Tikhonov smoothing matrix -----------------------------------%
switch order
    case 0 % 0th order Tikhonov
        if isa(n_grid, 'Grid'); Lpr0 = -speye(n_grid.Ne); % use Grid method (for partial grid support)
        else; Lpr0 = -speye(x_length);
        end
        
    case 1  % 1st order Tikhonov
        if isa(n_grid, 'Grid'); Lpr0 = n_grid.l1([], bc); % use Grid method (for partial grid support)
        else
            I1 = 0.5 .* speye(n, n);
            E1 = full(sparse(1:n-1, 2:n, 1, n, n));
            D1 = E1-I1;

            m = x_length / n;
            I2 = 0.5 .* speye(m, m);
            E2 = sparse(1:m-1, 2:m, 1, m, m);
            D2 = E2 - I2;

            Lpr0 = kron(I2, D1) + kron(D2, I1);

            Lpr0 = Lpr0 - spdiags(sum(Lpr0,2), 0, x_length, x_length);
            Lpr0(end,:) = [];
        end
        
    case 2  % 2nd order Tikhonov
        switch var{1}  % adds variants of 2nd order Tikhonov
            case 0  % standard Laplacian
                if isa(n_grid, 'Grid'); Lpr0 = n_grid.l2(bc);  % use Grid method (for partial grid support)
                else
                    I1 = 0.25 .* speye(n,n);
                    E1 = sparse(1:n-1,2:n,1,n,n);
                    D1 = E1 + E1' - I1;
        
                    m = x_length / n;
                    I2 = 0.25.*speye(m, m);
                    E2 = sparse(1:m-1, 2:m, 1, m, m);
                    D2 = E2+E2'-I2;
                    
                    Lpr1 = kron(I2, D1);
                    Lpr2 = kron(D2, I1);
                    Lpr0 = Lpr1 + Lpr2;
                    Lpr0 = Lpr0 - spdiags(sum(Lpr0,2), 0, x_length, x_length);
                end
        
            case 1  % difference in both dimensions
                if isa(n_grid, 'Grid')
                    n = n_grid.ne(1);
                    x_length = prod(n_grid.ne);  % for full gird
                end

                if length(var) == 1; var{2} = 1; end
                
                m = x_length / n;
                
                I1 = speye(n, n);
                D1 = -2 .* speye(m, m);
                D1 = spdiags(ones(m,2), -1, D1);
                D1 = spdiags(ones(m,2), 1, D1);
                
                I2 = speye(m, m);
                D2 = -2 .* speye(n, n);
                D2 = spdiags(ones(n,2), -1, D2);
                D2 = spdiags(ones(n,2), 1, D2);
        
                Lpr1 = var{2} .* kron(I2, D2);
                Lpr2 = kron(D1, I1);
                Lpr0 = Lpr1 - Lpr2;

                if isa(n_grid, 'PartialGrid')
                    Lpr0(n_grid.missing, :) = [];
                    Lpr0(:, n_grid.missing) = [];
                end
        end

    case 3  % 3rd order derivative
        if isa(n_grid, 'Grid')
            n = n_grid.ne(1);
            x_length = prod(n_grid.ne);  % for full gird
        end

        I1 = speye(n, n);
        m = x_length / n;
        
        D1 = sparse(m, m);
        D1 = spdiags(-1/2 .* ones(m,2), -2, D1);
        D1 = spdiags(ones(m,2), -1, D1);
        D1 = spdiags(-ones(m,2), 1, D1);
        D1 = spdiags(1/2 .* ones(m,2), 2, D1);
        % D1(1:2,:) = [];
        % D1(end-1:end,:) = [];
        
        I2 = speye(m, m);
        D2 = sparse(n, n);
        D2 = spdiags(-1/2 .* ones(n,2), -2, D2);
        D2 = spdiags(ones(n,2), -1, D2);
        D2 = spdiags(-ones(n,2), 1, D2);
        D2 = spdiags(1/2 .* ones(n,2), 2, D2);
        % D2(1:2,:) = [];
        % D2(end-1:end,:) = [];
        
        Lpr1 = kron(I2, D2);
        Lpr2 = kron(D1, I1);
        % Lpr0 = [Lpr1; Lpr2];
        Lpr0 = Lpr1 + Lpr2;

        if isa(n_grid, 'PartialGrid')
            Lpr0(n_grid.missing, :) = [];
            Lpr0(:, n_grid.missing) = [];
        end

    otherwise
        if and(order > 1, order < 2)  % 1st order Tikhonov with rotated matrix
            slope = (order - 1) * 10;
            if isa(n_grid,'Grid'); Lpr0 = n_grid.l1(slope, bc); % use Grid method (for partial grid support)
            else
                I1 = 1/2 .* speye(n, n);
                E1 = full(sparse(1:n-1, 2:n, 1, n, n));
                D1 = E1-I1;
        
                m = x_length / n;
                I2 = slope/2 .* speye(m, m);
                E2 = sparse(1:m-1, 2:m, 1, m, m);
                D2 = E2 - I2;
        
                Lpr0 = kron(I2, D1) + kron(D2, I1);
        
                Lpr0 = Lpr0 - spdiags(sum(Lpr0,2), 0, x_length, x_length);
                Lpr0(end,:) = [];
            end

        else
            disp('The specified order of Tikhonov is not available.');
            disp(' ');
        end
end

end


