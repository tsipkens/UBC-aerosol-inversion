
% GRID  Responsible for discretizing space as a grid and related operations.
%  
%  Grid is a class developed to discretize a parameter space (e.g., 
%  mass-mobility space). This is done using a simple rectangular grid that 
%  can have linear, logarithmic or custom spaced elements along the edges. 
%  Methods are designed to make it easier to deal with gridded data, 
%  allowing users to reshape vectorized data back to a 2D grid 
%  (`Grid.reshape` method) or vice versa. Other methods allow for plotting 
%  the 2D representation of vector data (`Grid.plot2d` method) or 
%  calculate the gradient of vector data (`Grid.grad` method).
%  
%  G = Grid(SPAN,NE) creates a grid with the domain specified by SPAN, a
%  2x2 array with [min(dim1),max(dim1); min(dim2),max(dim2)], and with the
%  number of elements in each dimension specified by NE, a 1x2 array. 
% 
%  G = Grid(EDGES) create a grid with edges specified by the entries of
%  EDGES, a 1x2 cell with the edges for the dim1 and dim2, respectively. 
%  Some functionality is limited if the edges are not uniform in log or 
%  linear space. 
% 
%  G = GRID(SPAN,NE,DISCRETE) adds an input to specify whether logarithmic
%  (default), specified using DISCRETE = 'log', or linear spacing, specified 
%  using DISCRETE = 'linear'. 
% 
%  AUTHOR: Timothy Sipkens, 2019-02-03
%  
%  ------------------------------------------------------------------------
%  
%  Instances of the Grid class can primarily be constructed in two ways. 
%  First, one can specify a `Grid.span` for the grid to cover in the 
%  parameter space. The span is specified using a 2 x 2 matrix, where the 
%  first row corresponds to the span for the first dimension of the 
%  parameter space (e.g., mass) and the second row corresponds to the span 
%  for the second dimension of the parameter space (e.g., mobility 
%  diameter). For example, if one wanted to logarithmically discretize 
%  mass space between 0.01 and 100 fg and mobility space between 10 and 
%  1000 nm, one could call:
%  
%  ```Matlab
%  span = [0.01,100; 10,1000]; % span of space to be covered
%  ne = [10,12]; % number of elements for each dimension
%  grid = Grid(span, ne, 'log'); % create instance of Grid
%  ```
%  
%  Second, one can supply a 1 x 2 cell array of edges, where the first 
%  entry is the center of the elements in the first dimension of parameter 
%  space and the second entry of the elements in the second dimension of 
%  parameter space. For example, to make a simple grid with elements at 
%  0.1 and 1 fg in mass space and 10, 200, and 1000 nm in mobility space, 
%  one would call:
%  
%  ```Matlab
%  edges = {[0.1,1], [10,200,1000]}; % cell array of edge vectors
%  grid = Grid(edges, [], 'log'); % create instance of Grid
%  ```
%  
%  Note that the number of elements is not required in this instance, as 
%  it is implied by the length of the vectors given in `edges`. The 
%  `'log'` (or equivalently `'logarithm'`) argument is still required to 
%  specify where nodes would be placed between the elements.
%  
%  Both the data, **b**, and two-dimensional size distribution, **x**, 
%  vectors can be defined with respect to an instance of this class. 
%  Generally, the data will only rely on the center of the elements on the 
%  grid (the width of the grid elements has little meaning for data). 
%  The vectors are arranged such that the first entry corresponds to the 
%  smallest size in both dimensions. The vector proceeds, first with 
%  increasing the first size dimension (e.g., for mass-mobility 
%  distributions this is mass by default) and then with increasing the 
%  second size dimension. Vectorizing the 2D gridded data can be done 
%  using the colon operand, i.e., `x(:)`, or using the `Grid.vectorize` 
%  method.
%  
%  For information on partial grids (where some elements are ignored), 
%  refer to `help PartialGrid` and `help Grid.partial`.

classdef Grid

properties
    discrete = {'log', 'log'};
                % type discretization to be applied to the edges
                % ('log'/'logarithmic' or 'linear')
    
    dim = 2;    % number of dimensions of mesh
    type = {};  % type of particle size on each dimension (e.g., mass, mobility)
    
    span = [];  % span covered by the grid in each dimension
                % Span applies to the center of the elements, i.e.
                % span(1,1) is the center of the first element for the
                % first dimension.
    
    ne = [];    % number of pixels/elements in each dimenion
    Ne = [];    % total number of pixels/elements, initially Ne = prod(ne)
                % reduced for partial grids
    
    edges = []; % cell of vectors containing edge points of pixel/element
                % centers, with one cell entry per dimension
    nodes = []; % contains position of nodes surrounding elements for each dimension
                % each cell has a vector of size (ne + 1).
    
    elements = [];  % contains position of pixel/element centers as a (ne x 2) vector
    nelements = []; % position of pixel/elements edges as a (ne x 4) vector
                    % [dim1_low,dim1_high,dim2_low,dim2_high]

    nodes_tr = [];  % transformed nodes to allow for log/lin basis
    nelements_tr = [];  % transformed nelements ...
                
    adj = [];   % adjacency matrix
end


methods
    %== GRID =========================================================%
    function obj = Grid(span_edges, ne, discrete)
    % GRID  Class constructor.
        
        %-- Parse inputs ---------------------------------------------%
        if nargin==0; return; end % return empty grid
        
        % If discrete is not specified, use logarithmic spacing.
        if ~exist('discrete','var'); discrete = []; end
        if isempty(discrete); discrete = 'log'; end
        if strcmp(discrete, 'logarithmic'); discrete = 'log'; end % allow for longhand
        if strcmp(discrete, 'linear'); discrete = 'lin'; end % allow for longhand
        if ischar(discrete); discrete = {discrete, discrete}; end % stretch for both dims
        
        if isa(span_edges,'cell')  % consider case where edges are given
            obj.edges = span_edges;
            obj.ne = [length(span_edges{1}),...
                length(span_edges{2})];
            obj.span = [min(span_edges{1}),max(span_edges{1});...
                min(span_edges{2}),max(span_edges{2})];

        else  % otherwise, consider case where span is given
            obj.span = span_edges;
            obj.ne = ne;
        end
        obj.Ne = prod(obj.ne);

        obj.discrete = discrete;

        obj = obj.mesh; % generates grid points
        obj = obj.adjacency; % get adjacency matrix
    end
    %=================================================================%
    
    

    %== MESH =========================================================%
    function obj = mesh(obj)
    % MESH  Responsible for generating a mesh represented by a series of elements.
    %  AUTHOR: Timothy Sipkens, 2019-02-03
    %  
    %  Currently setup to do simple linear or logarithmic spaced
    %  quadrilateral mesh.
    %
    %  obj.nodes contains the position of each of the nodes as
    %  a matrix, with a row for each node and a column for each
    %  dimension.
    
        %-- If required, generate edge discretization vectors --------%
        if isempty(obj.edges)
            for ii=1:obj.dim % loop through both dimensions
                if strcmp('lin',obj.discrete{ii})
                    obj.edges{ii} = linspace(obj.span(ii,1),obj.span(ii,2),obj.ne(ii));

                elseif strcmp('log',obj.discrete{ii})
                    obj.edges{ii} = logspace(...
                        log10(obj.span(ii,1)),log10(obj.span(ii,2)),obj.ne(ii));
                end
            end
            obj.Ne = prod(obj.ne);
        end

        %-- Generate nodes -------------------------------------------%
        for ii=1:obj.dim
            if strcmp(obj.discrete{ii},'log')
                r_m = exp((log(obj.edges{ii}(2:end))+...
                    log(obj.edges{ii}(1:(end-1))))./2); % mean of edges

                obj.nodes{ii} = [exp(2*log(obj.edges{ii}(1))-log(r_m(1))),...
                    r_m, exp(2*log(obj.edges{ii}(end))-log(r_m(end)))];
                
            elseif strcmp(obj.discrete{ii},'lin')
                r_m = (obj.edges{ii}(2:end)+...
                    obj.edges{ii}(1:(end-1)))./2; % mean of edges

                obj.nodes{ii} = [2*obj.edges{ii}(1)-r_m(1),...
                    r_m, 2*obj.edges{ii}(end)-r_m(end)];
            end
        end

        %-- Generate vectorized lists of elements --------------------%
        %   One column per dimension
        [vec1{1},vec1{2}] = ndgrid(obj.edges{1},obj.edges{2});
        obj.elements(:,1) = vec1{1}(:); % vectorize output
        obj.elements(:,2) = vec1{2}(:);
        
        [vec1{1},vec1{2}] = ndgrid(obj.nodes{1}(1:(end-1)),...
            obj.nodes{2}(1:(end-1)));
        [vec2{1},vec2{2}] = ndgrid(obj.nodes{1}(2:end),...
            obj.nodes{2}(2:end));
        obj.nelements = [vec1{1}(:),vec2{1}(:),vec1{2}(:),vec2{2}(:)];
    end
    %=================================================================%
    
    
    
    %== ADJACENCY ====================================================%
    function [obj, adj] = adjacency(obj, w)
    % ADJACENCY  Compute the adjacency matrix for the full grid using a four-point stencil.
    %  W is an optional that adds a weight for vertical pixels.
        
        if ~exist('w', 'var'); w = []; end
        if isempty(w); w = 1; end
        
        ind1 = ones(3 * prod(obj.ne), 1);
        ind2 = ones(3 * prod(obj.ne), 1);
        vec = zeros(3 * prod(obj.ne), 1);
        ll = 0;
        
        for jj=1:prod(obj.ne)
            if ~(mod(jj, obj.ne(1))==0) % up pixels
                ll = ll + 1;
                ind1(ll) = jj;
                ind2(ll) = jj + 1;
                vec(ll) = w;
            end
            
            if ~(mod(jj-1,obj.ne(1))==0) % down pixels
                ll = ll + 1;
                ind1(ll) = jj;
                ind2(ll) = jj - 1;
                vec(ll) = w;
            end
            
            if jj>obj.ne(1) % left pixels
                ll = ll + 1;
                ind1(ll) = jj;
                ind2(ll) = jj - obj.ne(1);
                vec(ll) = 1;
            end
            
            if jj <= (prod(obj.ne) - obj.ne(1)) % right pixels
                ll = ll + 1;
                ind1(ll) = jj;
                ind2(ll) = jj + obj.ne(1);
                vec(ll) = 1;
            end
        end
        
        adj = sparse(ind1, ind2, vec,...
            prod(obj.ne), prod(obj.ne));
        
        if isempty(obj.adj)  % for first run
            obj.adj = adj;
        end
    end
    %=================================================================%
    
    
    
    %== ADJACENCY8 ===================================================%
    function [adj] = adjacency8(obj)
    % ADJACENCY8  Compute the adjacency matrix for the full grid using an eight-point stencil.
        
        ind1 = ones(7*prod(obj.ne),1);
        ind2 = ones(7*prod(obj.ne),1);
        vec = zeros(7*prod(obj.ne),1);
        ll = 0;
        
        for jj=1:prod(obj.ne)
            if ~(mod(jj,obj.ne(1))==0) % up pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1;
                vec(ll) = 1;
            end
            
            if ~(mod(jj-1,obj.ne(1))==0) % down pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1;
                vec(ll) = 1;
            end
            
            if jj>obj.ne(1) % left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-obj.ne(1);
                vec(ll) = 1;
            end
            
            if jj<=(obj.Ne-obj.ne(1)) % right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj,obj.ne(1))==0),...
                    jj>obj.ne(1)) % up, left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1-obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj,obj.ne(1))==0),...
                    jj<=(obj.Ne-obj.ne(1))) % up, right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj+1+obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj-1,obj.ne(1))==0),...
                    jj>obj.ne(1)) % down, left pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1-obj.ne(1);
                vec(ll) = 1;
            end
            
            if and(~(mod(jj-1,obj.ne(1))==0),...
                    jj<=(obj.Ne-obj.ne(1))) % down, right pixels
                ll = ll+1;
                ind1(ll) = jj;
                ind2(ll) = jj-1+obj.ne(1);
                vec(ll) = 1;
            end
        end
        
        adj = sparse(ind1,ind2,vec,...
            prod(obj.ne),prod(obj.ne));
        
    end
    %=================================================================%
    
    
    
    %== GLOBAL_IDX ===================================================%
    function k = global_idx(obj, idx1, idx2)
    % GLOBAL_IDX  Convert 2D grid coordinate to a global coordinate in the grid.
    %  For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %  idx2 is the mobility index.
        
        k = idx1+(idx2-1)*obj.ne(1);
        
    end
    %=================================================================%
    
    
    
    %== TWO_IDX ======================================================%
    function [idx1,idx2] = two_idx(obj,k)
    % TWO_IDX  Convert global grid coordinate to 2D index on grid.
    %  For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %  idx2 is the mobility index.
    
        idx1 = mod(k,obj.ne(1));
        idx1(idx1==0) = obj.ne(1);
        
        idx2 = floor((k-1)./obj.ne(1))+1;
    end
    %=================================================================%
    
    
    
    %== L1 ===========================================================%
    function l1 = l1(obj, w, bc)
    % L1  Compute the first-order Tikhonov operator.
    %  W adds a weight used to reevaluate the adjacency matrix.
        
        if ~exist('bc', 'var'); bc = []; end
        if ~exist('w', 'var'); w = []; end
        if ~isempty(w)
            [~, adj_local] = obj.adjacency(w);  % re-evaluate adjacency with weight
        else
            adj_local = obj.adj;
        end
    
        l1 = -diag(sum(tril(adj_local))) + ...
            triu(adj_local);

        % Add unity on diagonal in final row for stability in square matrix
        % alternatively, this row can be deleted, however this causes
        % issues during Tikhonov inversion using this method
        l1(size(adj_local, 1), end) = -1;
        
        if bc == 0  % force zeros at BC. (update Jan. 2024)
            isedge = find(obj.elements(:,2) == obj.edges{2}(end));
            isedge = [isedge; find(obj.elements(:,1) == obj.edges{1}(end))];
            l1(isedge, isedge) = eye(length(isedge));  % replace entries with identity
        end
        
    end
    %=================================================================%
    
    

    %== L2 ===========================================================%
    function l2 = l2(obj, bc)
    % L2  Compute Laplacian
    %  Form is equiavalent to Tikhonov operator applying no slope at grid boundary.

        if ~exist('bc', 'var'); bc = 2; end

        l2 = -diag(sum(obj.adj)) + ...
            triu(obj.adj) + tril(obj.adj);
        
        if bc == 0  % force zeros at BC. (update Jan. 2024)
            isedge = find(obj.elements(:,2) == obj.edges{2}(end));
            isedge = [isedge; find(obj.elements(:,1) == obj.edges{1}(end))];
            isedge = [isedge; find(obj.elements(:,2) == obj.edges{2}(1))];
            isedge = [isedge; find(obj.elements(:,1) == obj.edges{1}(1))];
            l2(isedge, isedge) = eye(length(isedge));  % replace entries with identity
        end
    end
    %=================================================================%
    
    
    
    %== GRAD =========================================================%
    function out = grad(obj,x)
    % GRAD  Calculates the gradient in x. 
    %  Uses simple first order differences.

        [~,dr1,dr2] = obj.dr;

        [grad2,grad1] = gradient(reshape(x,obj.ne));
        grad1 = grad1.*dr1;
        grad2 = grad2.*dr2;

        out = [grad1(:),grad2(:)];

    end
    %=================================================================%
    
    
    
    %== PROJECT ======================================================%
    function x = project(obj,grid_old,x)
    % PROJECT  Project x onto current grid. 
    %  Uses simple linear interpolation for this purpose. The parameter 
    %  GRID_OLD contains the original grid for the input data X.
        
        if isa(grid_old, 'PartialGrid')  % added processing for partial grids
            x = grid_old.partial2full(x);
        end
        
        n1 = length(grid_old.edges{1});
        n2 = length(grid_old.edges{2});
        x = reshape(x,[n1,n2]);
        
        [edges1,edges2] = ndgrid(obj.edges{1},obj.edges{2});
        [edges_old1,edges_old2] = ndgrid(grid_old.edges{1},grid_old.edges{2});
        F = griddedInterpolant(edges_old1,edges_old2,x,'linear','linear');
        
        x = F(edges1,edges2);
        x = x(:);
        
    end
    %=================================================================%



    %== TRANSFORM ====================================================%
    function B = transform(obj, grid_old)
    % TRANSFORM  Function to transform kernel functions. Output is a matrix to 
    %   be multiplied by the original kernel, A.
        
        for ii=1:obj.dim % loop over both dimensions
            dr_inv = 1./(grid_old.nodes{ii}(2:end)-grid_old.nodes{ii}(1:(end-1)));

            t0 = min(1,max(0,...
                bsxfun(@times,...
                grid_old.nodes{ii}(2:end)-obj.nodes{ii}(1:(end-1))',dr_inv)...
                )); % upper matrix of overlap
            t1 = min(1,max(0,...
                bsxfun(@times,...
                obj.nodes{ii}(2:end)'-grid_old.nodes{ii}(1:(end-1)),dr_inv)...
                )); % lower matrix of overlap
            t2{ii} = sparse(min(t0,t1));
        end % end loop over both dimensions

        ind_old_tot = 1:prod(grid_old.ne);
        ind1_old = mod(ind_old_tot-1,grid_old.ne(1))+1;
        ind2_old = ceil(ind_old_tot./grid_old.ne(1));
        
        ind_tot = 1:prod(obj.ne);
        ind1 = mod(ind_tot-1,obj.ne(1))+1;
        ind2 = ceil(ind_tot./obj.ne(1));
        
        B = (t2{1}(ind1,ind1_old).*t2{2}(ind2,ind2_old))';
    end
    %=================================================================%



    %== DR ===========================================================%
    function [dr,dr1,dr2] = dr(obj)
    % DR  Calculates the differential area of the elements in the grid.
        
        dr_0 = cell(obj.dim, 1);
        for ii=1:obj.dim
            if strcmp(obj.discrete{ii}, 'log')
                dr_0{ii} = log10(obj.nodes{ii}(2:end)) - ...
                    log10(obj.nodes{ii}(1:(end-1)));
            
            elseif strcmp(obj.discrete{ii},'lin')
                dr_0{ii} = obj.nodes{ii}(2:end) - ...
                    obj.nodes{ii}(1:(end-1));
            end
        end
        
        [dr1,dr2] = ndgrid(dr_0{1}, dr_0{2});
        dr1 = abs(dr1); % in case edges vector is reversed
        dr2 = abs(dr2);
        dr = dr1(:) .* dr2(:);
        
    end
    %=================================================================%


    
    %== MARGINALIZE ==================================================%
    function [marg, tot] = marginalize(obj, x, dim)
    % MARGINALIZE  Marginalizes over the grid in each dimension.
    %   Uses Euler's method to integrate over domain.
        
        x = obj.reshape(x);
        
        [dr,dr1,dr2] = obj.dr; % generate differential area of elements
        dr = obj.reshape(dr); % fills out dr for partial grids
        
        tot = sum(x(:).*dr(:)); % integrated total
        
        t0 = dr2; % added processing for partial elements
        dr2 = dr./dr1;
        dr1 = dr./t0;
        
        marg{1} = sum(dr2 .* x,2); % integrate over diameter
        marg{2} = sum(dr1 .* x,1); % integrate over mass
        
        % If dim input, output marginalized distribution for 
        % specific dimension. 
        if exist('dim', 'var')
            marg = marg{dim};
        end
        
    end
    %=================================================================%
    
    
    
    %== MARGINALIZE_OP ===============================================%
    function [C1, dr0] = marginalize_op(obj,dim)
    % MARGINALIZE_OP  A marginalizing operator, C1, to act on 2D distributions.
    %  AUTHOR: Timothy Sipkens, Arash Naseri, 2020-03-09
        
        if ~exist('dim','var'); dim = []; end
        if isempty(dim); dim = 1; end
            
        switch dim % determine which dimension to sum over
            case 2 % integrate in column direction
                C1 = kron(speye(obj.ne(2)),ones(1,obj.ne(1)));

            case 1 % integrate in row direction
                C1 = repmat(speye(obj.ne(1),obj.ne(1)),[1,obj.ne(2)]);
        end
        
        dr0 = ones(obj.Ne,1);
        
    end
    
    
    
    %== RESHAPE ======================================================%
    function x = reshape(obj, x)
    % RESHAPE  A simple utility to reshape a vector based on the grid.
    %  X = Grid.reshape(X) reshapes X based on the structure of Grid.
        x = reshape(x, obj.ne);
    end
    %=================================================================%



    %== VECTORIZE ====================================================%
    function [x,t1,t2] = vectorize(obj,x)
    % VECTORIZE  A simple function to vectorize 2D gridded data.
    % 
    %  OUTPUTS:
    %   x	Vectorized data
    %   t1  Vectorized element centers for first dimension
    %   t2  Vectorized element centers for second dimension
    
        if exist('x','var'); x = x(:); else; x = []; end
        
        if nargout>1; t1 = obj.elements(:,1); end
        if nargout>2; t2 = obj.elements(:,2); end
    end
    %=================================================================%



    %== RAY_SUM ======================================================%
    function [C,rmin,rmax] = ray_sum(obj,logr0,slope,f_bar)
    % RAY_SUM  Perfrom a ray sum for a given ray and the current grid.
    %   Currently assumes uniform, logarithmic grid 
    %   and can accommodate partial grids.
    %   BASIS:  Code from Samuel Grauer
    %   AUTHOR: Timothy Sipkens, 2019-07-14
    % 
    %  INPUTS:
    %   logr0   A single point on the line in log-log space, r0 = log10([dim1,dim2])
    %	slope   Slope of the line
    %   f_bar   Flag for progress bar
    %  
    %  OUTPUTS:
    %   C       Ray-sum matrix
        
        if ~exist('f_bar','var'); f_bar = []; end
        if isempty(f_bar); f_bar = 0; end
        
        %-- Preallocate arrays ---------------------------------------%
        m = size(slope,1);
        C = spalloc(m,obj.Ne,ceil(0.1*m*obj.Ne)); % assume 10% full
        
        
        %-- Compute ray-sum matrix -----------------------------------%
        if f_bar; tools.textbar(); end
        for ii=1:m % loop over multiple rays

            %-- Ray vector -------------%
            dv = [1,slope]; % convert slope to step vector along line
            dv = dv/norm(dv);
            dv(dv == 0) = 1e-10; % for stability during division
            
            
            %-- Line intersections -----%
            %   Use parametric representation of the line and finds two
            %   intersections for each element.
            %   Note: Assumes a logarithmic grid.
            tmin = (log10(obj.nelements(:,[3,1]))-...
                logr0)./dv; % minimum of element
            tmax = (log10(obj.nelements(:,[4,2]))-...
                logr0)./dv; % maximum of element
            
            
            %-- Corrections ------------%
            %   Decide which points would correspond to transecting the
            %   pixel.
            tmin = max(tmin,[],2);
            tmax = min(tmax,[],2);
            
            
            %-- Convert back to [x,y] --%
            rmin = logr0 + tmin.*dv; % location of intersect with min. of pixel
            rmax = logr0 + tmax.*dv; % location of intersect with max. of pixel
            
            rmin = min(rmin, log10(obj.nelements(:, [4,2])));
            rmax = max(rmax, log10(obj.nelements(:, [3,1])));
            
            chord = sqrt(sum((rmax - rmin).^2,2)); % chord length
            chord(chord<1e-15) = 0; % truncate small values
            
            
            %-- Ray-sum matrix ---------%
            [~,jj,a] = find(chord');
            if ~isempty(a)
                C(ii,:) = sparse(1, jj, a, 1, obj.Ne, ceil(obj.Ne));
            end
            if f_bar, tools.textbar(ii/m); end
            
            
            %-- Modify rmin and rmax for output -----%
            rmin = fliplr(rmin);
            rmax = fliplr(rmax);
            
        end % end loop over multiple rays

    end
    %=================================================================%
    
    
    
    %== CLOSEST_IDX =================================================%
    function [k,idx_2d] = closest_idx(obj,r0)
    % CLOSEST_IDX  Returns the pixel in which r0 is located.
    %   This function uses vector operations to find multiple points.
    % 
    %  INPUTS:
    %   r0      Coordinates in grid space, r0 = [dim1,dim2]
    %           Can form N x 2 vector, where N is the number of points
    %           to be found.
    %  OUTPUTS:
    %   k       Global index on the grid, incorporating missing pixels
    %   idx_2d  Pair of indices of pixel location
        
        idx_2d = zeros(size(r0,1),2); % pre-allocate
        
        for ii=1:obj.dim
            idx_bool = and(...
                r0(:,ii)>=obj.nodes{ii}(1:(end-1)),...
                r0(:,ii)<obj.nodes{ii}(2:end));
                % if on border, place in higher pixel
            
            [~,idx_2d(:,ii)] = max(idx_bool,[],2);
            
            idx_2d(~any(idx_bool,2),ii) = NaN;
                % if point is outside of grid, return NaN
        end
        
        f_nan = or(isnan(idx_2d(:,1)),isnan(idx_2d(:,2)));
        k = NaN(size(r0,1),1);
        k(~f_nan) = obj.global_idx(idx_2d(~f_nan,1),idx_2d(~f_nan,2));
    end
    %=================================================================%

    
    
%=====================================================================%
%-- VISUALIZATION METHODS --------------------------------------------%
%=====================================================================%
    
    %== PLOT2D =======================================================%
    function [h,x] = plot2d(obj, x, f_tr, f_contf)
    % PLOT2D  Plots x as a 2D function on the grid.
    %   AUTHOR: Timothy Sipkens, 2018-11-21
        
        if ~exist('f_contf','var'); f_contf = []; end % set empty contourf flag
        if isempty(f_contf); f_contf = 0; end % set contourf flag to false
        
        if ~exist('f_tr','var'); f_tr = []; end
        if isempty(f_tr); f_tr = 0; end % set transpose flag to false
        
        cla; % clear existing axis
        
        %-- Issue warning if grid edges are not be uniform -----------%
        %   The imagesc function used here does not conserve proportions.
        [~, dr1, dr2] = obj.dr;  % used to give warning below
        dr0 = dr1(:) .* dr2(:);
        if ~all(abs(dr0(2:end) - dr0(1)) < 1e-10)
            warning(['The plot2d method does not necessarily display ',...
                'correct proportions for non-uniform grids. ', ...
                'Use a regularily-spaced grid for best results.']);
        end
        
        x = obj.reshape(x);  % reshape

        % Default for order of indices. Plus transpose if flagged.
        idx1 = 1;  idx2 = 2;
        if f_tr; idx1 = 2; idx2 = 1; x = x'; end
        
        %-- Plot -------------------------------%
        if f_contf==0 % plot as image
            imagesc(obj.edges{idx2}, obj.edges{idx1},x);
            set(gca,'YDir','normal');
        else % plot as contourf
            contourf(obj.edges{idx2},obj.edges{idx1},x,35,'EdgeColor','none');
        end
        
        %-- Adjust tick marks for log scale ----%
        %   'logarithmic' included for backward compatibility.
        if strcmp(obj.discrete{2}, 'log')
            set(gca, 'XScale', 'log');
        end
        if strcmp(obj.discrete{1}, 'log')
            set(gca, 'YScale', 'log');
        end
        
        xlim(obj.span(idx2,:));  % set plot limits bsed on grid limits
        ylim(obj.span(idx1,:));

        if ~isempty(obj.type)
            xlabel(obj.type{idx2});
            ylabel(obj.type{idx1});
        end
        
        % Grey labels and axes to allow viz against dark and light bgs.
        set(gca, 'XColor', [0.5, 0.5, 0.5], ...
            'YColor', [0.5, 0.5, 0.5], ...
            'linewidth', 0.75);

        axis square;
        
        if nargout>0; h = gca; end
        
    end
    %=================================================================%



    %== PLOT2D_MARG ==================================================%
    function [h,x_m] = plot2d_marg(obj,x,obj_t,x_t,f_contf)
    % PLOT2D_MARG  Plots x as a 2D function on the grid, with marginalized distributions.
    %   AUTHOR: Timothy Sipkens, 2018-11-21
        
        if ~exist('f_contf','var'); f_contf = []; end % set empty contourf flag
        if ~exist('x_t','var'); x_t = []; end
        
        subplot(4,4,[5,15]);
        obj.plot2d(x,f_contf);
        
        x_m = obj.marginalize(x);
        
        
        %-- Plot marginal distribution (dim 2) -----------------------%
        subplot(4,4,[1,3]);
        marg_dim = 2;
        stairs(obj.nodes{marg_dim},...
            [x_m{marg_dim},0],'k');
        xlim([min(obj.edges{marg_dim}),max(obj.edges{marg_dim})]);
        set(gca,'XScale','log');
        
        %-- Also plot marginal of the true distribution --------------%
        if ~isempty(x_t)
            x_m_t = obj_t.marginalize(x_t);
            
            hold on;
            plot(obj_t.nodes{marg_dim},...
                [x_m_t{marg_dim},0],'color',[0.6,0.6,0.6]);
            hold off;
        end
        
        
        %-- Plot marginal distribution (dim 1) -----------------------%
        subplot(4,4,[8,16]);
        marg_dim = 1;
        stairs([0;x_m{marg_dim}],...
            obj.nodes{marg_dim},'k');
        ylim([min(obj.edges{marg_dim}),max(obj.edges{marg_dim})]);
        set(gca,'YScale','log');
        
        %-- Also plot marginal of the true distribution --------------%
        if ~isempty(x_t)
            hold on;
            plot([0;x_m_t{marg_dim}],...
                obj_t.nodes{marg_dim},'color',[0.6,0.6,0.6]);
            hold off;
        end
        
        subplot(4,4,[5,15]);
        if nargout>0; h = gca; end
        
    end
    %=================================================================%



    %== PLOT2D_SWEEP =================================================%
    function [h,x] = plot2d_sweep(grid,x,cm,dim)
    % PLOT2D_SWEEP  Plot data in slices, sweeping through the provided colormap.
    %   AUTHOR: Timothy Sipkens, 2019-11-28
        
        if ~exist('dim','var'); dim = []; end
        if isempty(dim); dim = 1; end
            % dimension to sweep through
            % e.g. sweep through mass setpoints on standard grid, dim = 1
        
        if ~exist('cm','var'); cm = []; end
        if isempty(cm); cm = colormap('gray'); end
        
        
        dim2 = setdiff([1,2],dim); % other dimension, dimension to plot
        
        addpath('cmap'); % load cmap package to use `sweep_cmap(...)`
        if isfile('cmap/cmap_sweep.m'); cmap_sweep(grid.ne(dim), cm); % set color order to sweep through colormap
        else; warning('The `cmap` package missing.'); end % if package is missing
        
        x_rs = reshape(x, grid.ne);
        if dim==2; x_rs = x_rs'; end
        
        % h = semilogx(grid.edges{dim2},x_rs,...
        %     'o-','MarkerSize',2.5,'MarkerFaceColor',[1,1,1]);
        h = semilogx(grid.edges{dim2}, x_rs, '-');

        if nargout==0; clear h; end
    
    end
    %=================================================================%
    
    
    
    %== PLOT_MARGINAL ================================================%
    function [] = plot_marginal(obj,x,dim,x0)
    % PLOT_MARGINAL  Plot marginal distributions.
    %   AUTHOR:	Timothy Sipkens, 2019-07-17
    %   NOTE: 'x' can be a cell array containing multiple x vectors

        %-- Parse inputs ---------------------------------------------%
        if ~iscell(x); x = {x}; end
            % if input is not cell, covert it to one

        if ~exist('dim','var'); dim = []; end
        if ~exist('x0','var'); x0 = []; end

        if isempty(dim); dim = 1; end
        if isempty(x0); x0 = x{1}; end

        x0_m = obj.marginalize(x0); % reference case


        %-- Plot entries in x ----------------------------------------%
        for ii=1:length(x) % plot other provided x
            x_m = obj.marginalize(x{ii});

            %-- Plot difference --%
            subplot(3,1,1);
            if ~isempty(findall(gca,'type','line')); hold on; end
                % if not the first line in the plot, hold on
            semilogx(obj.edges{dim},x_m{dim}-x0_m{dim});
            hold off;

            %-- Plot marginal distribution --%
            subplot(3,1,2:3);
            if ~isempty(findall(gca,'type','line')); hold on; end
                % if not the first line in the plot, hold on
            semilogx(obj.edges{dim}, x_m{dim});
            hold off;
        end


        %-- Set axes limits ------------------------------------------%
        subplot(3,1,2:3);
        xlim([min(obj.edges{dim}), max(obj.edges{dim})]);

        subplot(3,1,1);
        xlim([min(obj.edges{dim}), max(obj.edges{dim})]);
    end
    %=================================================================%



    %== PLOT_CONDITIONAL =============================================%
    function [] = plot_conditional(obj, x, dim, ind, x0)
    % PLOT_CONDITIONAL  Plot conditional distributions
    %  AUTHOR:	Timothy Sipkens, 2019-07-17
    %  NOTE: 'x' can be a cell array containing multiple x vectors

        %-- Parse inputs ---------------------------------------------%
        if ~iscell(x); x = {x}; end
            % if input is not cell, covert it to one

        if ~exist('dim','var'); dim = []; end
        if ~exist('ind','var'); ind = []; end
        if ~exist('x0','var'); x0 = []; end

        if isempty(dim); dim = 1; end
        if isempty(ind); ind = round(obj.ne(dim)/2); end
        if isempty(x0); x0 = x{1}; end

        x0 = reshape(x0,obj.ne); % reference case
        if dim==1; x0_c = x0(:,ind); end
        if dim==2; x0_c = x0(ind,:); end


        %-- Plot entries in x ----------------------------------------%
        for ii=1:length(x) % plot other provided x
            x_c = reshape(x{ii},obj.ne);
            if dim==1; x_c = x_c(:,ind); end
            if dim==2; x_c = x_c(ind,:); end

            %-- Plot difference --%
            subplot(3,1,1);
            if ~isempty(findall(gca,'type','line')); hold on; end
                % if not the first line in the plot, hold on
            semilogx(obj.edges{dim},x_c-x0_c);
            hold off;

            %-- Plot marginal distribution --%
            subplot(3,1,2:3);
            if ~isempty(findall(gca,'type','line')); hold on; end
                % if not the first line in the plot, hold on
            semilogx(obj.edges{dim},x_c);
            hold off;
        end


        %-- Set axes limits ------------------------------------------%
        subplot(3,1,2:3);
        xlim([min(obj.edges{dim}),max(obj.edges{dim})]);

        subplot(3,1,1);
        xlim([min(obj.edges{dim}),max(obj.edges{dim})]);
    end
    %=================================================================%
    
    
    
    %== PLOT2D_SCATTER ===============================================%
    function [] = plot2d_scatter(obj, x, cm)
    % PLOT2D_SCATTER  Wrapper for tools.plot2d_scatter.
    %  AUTHOR: Timothy Sipkens, 2020-11-05
    %  NOTE: 'x' can be a cell array containing multiple x vectors
        
        if ~exist('cm', 'var'); cm = []; end

        %-- Parse inputs ---------------------------------------------%
        [e2, e1] = meshgrid(obj.edges{2}, obj.edges{1});
        tools.plot2d_scatter(e1(:), e2(:), x, cm);
    end
    %=================================================================%
    
    
    
    %== PLOT2D_SCATTER ===============================================%
    function [] = surf(obj, x, cm)
    % SURF  Wrapper for tools.plot2d_scatter.
    %  AUTHOR: Timothy Sipkens, 2021-03-30
        
        if ~exist('cm', 'var'); cm = gray(255); end

        %-- Parse inputs ---------------------------------------------%
        surf(obj.edges{2}, obj.edges{1}, obj.reshape(x));
        colormap(cm);
        shading interp;
        
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
        
        set(gca, 'View', [45, 70]);  % set axis view
        
        xlim([min(obj.edges{2}),max(obj.edges{2})]);
        ylim([min(obj.edges{1}),max(obj.edges{1})]);
    end
    %=================================================================%
    
    
    
%=====================================================================%
%-- SUPPORT FOR PARTIAL GRIDS ----------------------------------------%
%=====================================================================%
    
    %== PARTIAL ======================================================%
    function partialgrid = partial(obj, varargin)
    % PARTIAL  Convert grid to a partial grid, removing elements above or 
    %  below a line. 
    % 
    %  G = Grid.partial(~, R0, SLOPE0) removes elements above the line that goes
    %  through the point R0 and having a slope of SLOPE0. For logarithmic
    %  grids, lines correspond to exponential curves, R0 are given as 
    %  log10(...) quantities, and slopes correspond to the exponent. For 
    %  example, Grid.partial([0,2],3) removes all grid elements above the
    %  exponential curve that passes through [1,100] and having an exponent
    %  of 3 (e.g., curves that increase volumetrically). 
    % 
    %  G = Grid.partial(~, R0, SLOPE0, R1, SLOPE1) adds a second set of 
    %  arguments analogous to above but remove points below a given line 
    %  (instead of above). 
        
        partialgrid = PartialGrid(obj.edges, [], obj.discrete, varargin{:});
        partialgrid.type = obj.type;
        
    end
    %=================================================================%

end

end
