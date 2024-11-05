
% PARTIALGRID  A subclass of Grid handling partially truncated grids.
%  
%  G = PartialGrid(...,R0,SLOPE0) takes the same inputs as the Grid
%  class constructor, but adds extra inputs to form a partial grid. 
%  This call removes elements above the line that goes through the 
%  point R0 and having a slope of SLOPE0. For logarithmic
%  grids, lines correspond to exponential curves, R0 are given as 
%  log10(...) quantities, and slopes correspond to the exponent. For 
%  example, Grid.partial([0,2],3) removes all grid elements above the
%  exponential curve that passes through [1,100] and having an exponent
%  of 3 (e.g., curves that increase volumetrically). 
%  
%  G = PartialGrid(...,R0,SLOPE0,R1,SLOPE1) adds a second set of 
%  arguments analogous to above but remove points below a given line 
%  (instead of above). 
%  
%  AUTHOR: Timothy Sipkens, 2021-03-29
%  
%  ------------------------------------------------------------------------
%  
%  Partial grids can be created two ways. First, if one has already 
%  created an instance of the **Grid** class, one can generate a 
%  **PartialGrid** using the `Grid.partial(...)` method, which cuts the 
%  grid about a line specified by a y-intercept and slope. Second, one can 
%  call the `PartialGrid(...)` constructor method directly, which takes 
%  that same first three arguments as the **Grid** class and appends the 
%  arguments from the `Grid.partial(...)` method. In either case, all of 
%  the grid elements with a center above this line is removed from the 
%  grid. For example, if one wanted to create a grid where all of the 
%  points above the 1-1 line should be removed (as is relevant for PMA-SP2 
%  inversion), one can call
%  
%  ```Matlab
%  grid = PartialGrid(0, 1);
%  ```
%  
%  which replaces the current **Grid** with a **PartialGrid** with the 
%  same span. Equivalently, assuming one is using the **Grid** class, one 
%  can generate this **PartialGrid** using
%  
%  ```Matlab
%  span = [0.01,100; 10,1000]; % span of space to be covered
%  ne = [10,12]; % number of elements for each dimension
%  grid = PartialGrid(span, ne, 'log', 0, 1);
%  ```
%  
%  For partial grids:
%  
%  1. `Grid.missing` contains a list of the global indices for the missing pixels;
%  2. `Grid.cut` contains the y-intercept and slope of the line used in 
%     making the cut (which is used in marginalization to only partially 
%     weight points that are cut off);
%  3. `Grid.elements` is updated to only list the center of the elements 
%     that are left in the grid (i.e., have indices that are not listed in 
%     `Grid.ismissing`);
%  4. `Grid.nelements` similarly only lists the edges of pixels that remain 
%     in the partial grid; and
%  5. `Grid.adj` is updated to only list the elements/pixels that are 
%     adjacent on the partial grid.
%  
%  The `PartialGrid.edges` and `PartialGrid.nodes` are inherited from the 
%  parent **Grid** and refer to the underlying grid structure prior to any 
%  of the grid points being removed.
%  
%  Methods specific to partial grids include:
%  1. `PartialGrid.partial2full(x)` will convert a partial x (resolved 
%     with only the remaining points on the grid) to one resolved on the 
%     original grid, filling the missing points with zeros. I use this to 
%     plot and marginalize the distributions (i.e., plots will have zeros 
%     in all of the missing elements/pixels) and
%  2. `PartialGrid.full2partial(x)` will do the reverse of above, 
%     converting a full x (with all of the points from the original grid) 
%     to one resolved on the partial grid.
%  
%  Most other Grid methods will operate on partial grids. For example, 
%  the `PartialGrid.marginalization(...)` method will only sum the pixels 
%  that are still on the grid (exploiting the fact that the other pixels 
%  are necessarily zero), including accounting for the fact that pixels on 
%  the diagonal are only half pixels (i.e., they are triangles rather than 
%  rectangular elements). The `PartialGrid.plot2d(...)` method, equally, 
%  will plot zeros in the upper left triangle. The `PartialGrid.l1(...)` 
%  method will generate the first-order Tikhonov matrix, 
%  **L**<sub>pr</sub>, for the modified grid (only considering pixels that 
%  are adjacent on the partial grid). The list goes on, and we refer the 
%  reader to the class definition below to view the modified methods. 

classdef PartialGrid < Grid

properties
    
    missing = [];  % global indices of missing grid points for partial grids
    cut = [];      % [y-intercept,slope] used to cut partial grid
    
end 


methods
        
    %== PARTIALGRID ======================================================%
    function obj = PartialGrid(span_edges, ne, discrete, r0, slope0, r1, slope1)
    % PARTIALGRID  Class cosntructor.
        
        % Call Grid constructor.
        obj = obj@Grid(span_edges, ne, discrete);

        %-- Parse inputs ----------------------------%
        if ~exist('slope0','var'); slope0 = []; end
        if isempty(slope0); slope0 = 0; end
        
        if ~exist('r0','var'); r0 = []; end
        if length(r0)==1; b0 = r0; end % if scalar, use as y-intercept
        if length(r0)==2; b0 = r0(1) - slope0 * r0(2); end % if coordinates, find y-intercept
        if isempty(r0); b0 = inf; end % if not specified, use b = inf
        
        %-- For bottom triangle --%
        if ~exist('slope1','var'); slope1 = []; end
        if isempty(slope1); slope1 = 0; end

        if ~exist('r1','var'); r1 = -inf; end
        if length(r1)==1; b1 = r1; end % if scalar, use as y-intercept
        if length(r1)==2; b1 = r1(1) - slope1 * r1(2); end % if coordinates, find y-intercept
        if isempty(r1); b1 = -inf; end % if not specified, use b = -inf
        %--------------------------------------------%

        %-- Cut upper triangle ---------------------%
        if strcmp(obj.discrete,'log')
            tup = log10(obj.nelements(:,[1,4]));
        else
            tup = obj.nelements(:,[1,4]);
        end
        tup = tup + abs(1e-3 .* mean(tup(2:end,:) - tup(1:(end-1),:))) .* [1,0];
            % avoids minimially overlapping elements
        
        f_missing = tup(:,1) > (tup(:,2) .* slope0 + b0);
        t1 = 1:prod(obj.ne);
        obj.cut = [b0, slope0];


        %-- Consider cutting lower triangle --------%
        if strcmp(obj.discrete,'log')
            tlow = log10(obj.nelements(:,[2,3]));
        else
            tlow = obj.nelements(:,[2,3]);
        end
        tlow = tlow+abs(1e-3.*mean(tlow(2:end,:)-tlow(1:(end-1),:))).*[-1,0];
            % avoids minimially overlapping elements

        if ~isinf(r1)
            f_missing1 = tlow(:,1)<(tlow(:,2).*slope1+b1);
            f_missing = or(f_missing,f_missing1);
            obj.cut = [obj.cut,b1,slope1];
        end


        %-- Update grid properties -----------------%
        obj.missing = t1(f_missing);

        obj.elements = obj.elements(~f_missing,:);
        obj.nelements = obj.nelements(~f_missing,:);
        obj.Ne = size(obj.elements,1);
        obj = obj.adjacency;

    end
    
    
    
    %== ADJACENCY ====================================================%
    function [obj, adj, isedge] = adjacency(obj, w)
    % ADJACENCY  Compute the adjacency matrix for the full grid using a four-point stencil.
    %  W is an optional that adds a weight for vertical pixels.
        
        if ~exist('w', 'var'); w = []; end
        
        [~, adj] = adjacency@Grid(obj, w);

        % Flags whether element is next to new edge.
        adju = triu(adj, 2);
        isedge = any(adju(:, obj.missing), 2);
        isedge(obj.missing) = [];
        
        adj(obj.missing, :) = [];
        adj(:, obj.missing) = [];
        
        obj.adj = adj;
        
    end
    
    
    
    %== GLOBAL_IDX ===================================================%
    function k = global_idx(obj, idx1, idx2)
    % GLOBAL_IDX  Convert 2D grid coordinate to a global coordinate in the grid.
    %  For mass-mobiltiy grids, for example, idx1 is the mass index and 
    %  idx2 is the mobility index.
        
        k = global_idx@Grid(obj, idx1, idx2);  % use Grid function
        
        %-- Update for partial grid. --%
        b0 = ismember((1:prod(obj.ne))',obj.missing);
        t0 = cumsum(b0); % number of entries missing prior to current index

        [~,i0] = intersect(k,obj.missing); % find indices in k that are missing

        k = k-t0(k); % reduce indices based on missing indices
        k(i0) = NaN; % return NaN for indices that are missing
    end
    %=================================================================%



    %== L1 ===========================================================%
    function [l1] = l1(obj, w, bc)
    % L1  Bridge to compute the first-order Tikhonov operator.
    %  Handles zeroing BC for particle grid boundaries.
        
        if ~exist('w', 'var'); w = []; end
        if ~exist('bc', 'var'); bc = []; end
        
        l1 = l1@Grid(obj, w);  % use inherited function
        
        % Force zeros at BC.
        % Requires determining which pixels are adjacent to remove points.
        if bc == 0
            [~, ~, isedge] = obj.adjacency(w);  % get missing elements from adj
            isedge = find(isedge);
            l1(isedge, isedge) = eye(length(isedge));
        end
        
    end
    %=================================================================%%
    
    
    
    %== PROJECT ======================================================%
    function x = project(obj, grid_old, x)
    % PROJECT  Project x onto current grid. 
    %  Uses simple linear interpolation for this purpose. The parameter 
    %  GRID_OLD contains the original grid for the input data X.
        
        x = project@Grid(obj, grid_old, x);
        
        x = obj.full2partial(x);
        
    end
    %=================================================================%
    
    
    
    %== TRANSFORM ====================================================%
    function B = transform(obj, grid_old)
    % TRANFORM  Function to transform kernel functions. 
    %  Output is a matrix to be multiplied by the original kernel, A.
        
        B = transform@Grid(obj, grid_old);
        
        B(:,obj.missing) = [];
        if isa(grid_old, 'PartialGrid')
            B(grid_old.missing,:) = [];
        end
        %-------------------------------------------------------------%
    end
    %=================================================================%
    
    
    
    %== DR ===========================================================%
    function [dr, dr1, dr2] = dr(obj)
    % DR  Calculates the differential area of the elements in the grid.
        
        [dr, dr1, dr2] = dr@Grid(obj);
        
        %-- Added processing for partial grids -----------------------%
        dr = obj.full2partial(dr); % used if lower cut is employed
        
        %{
        % Archived code that could calculate remaining area after cut. 
        % When commented, full area of element if maintained if any part 
        % of the element remains. 
        [~,rmin,rmax] = obj.ray_sum([0, obj.cut(1)], obj.cut(2), 0);
        t0 = (rmin(:,1) - log10(obj.nelements(:,1))) .* ...
            (log10(obj.nelements(:,4)) - log10(obj.nelements(:,3)));
                % lower rectangle
        t1 = (log10(obj.nelements(:,4)) - rmax(:,2)).*...
            (log10(obj.nelements(:,2)) - rmin(:,1));
                % right rectangle
        t2 = 1/2 .* (rmax(:,1) - rmin(:,1)) .* ...
            (rmax(:,2) - rmin(:,2));
                % upper, left triangle
        dr = t0+t1+t2;

        if length(obj.cut)==4 % consider cutting lower triangle
            [~, rmin, rmax] = obj.ray_sum([0, obj.cut(3)], obj.cut(4), 0);
            t0 = (log10(obj.nelements(:,2)) - rmax(:,1)) .* ...
                (log10(obj.nelements(:,4)) - log10(obj.nelements(:,3)));
                    % upper rectangle
            t1 = (rmin(:,2) - log10(obj.nelements(:,3))) .* ...
                (rmax(:,1) - log10(obj.nelements(:,1)));
                    % left rectangle
            t2 = 1/2 .* (rmax(:,1) - rmin(:,1)).*...
                (rmax(:,2) - rmin(:,2));
                    % lower, right triangle
            dr = dr .* (t0+t1+t2) ./ dr0; % accounts for element that are discected twice
        end
        %}
        %-------------------------------------------------------------%
    end
    %=================================================================%
    
    
    
    %== MARGINALIZE_OP ===============================================%
    function [C1,dr0] = marginalize_op(obj, dim)
    % MARGINALIZE_OP  A marginalizing operator, C1, to act on 2D distributions.
    %  AUTHORS: Timothy Sipkens, Arash Naseri, 2020-03-09
        
        C1 = marginalize_op@Grid(obj, dim);
        
        C1(:,obj.missing) = []; % remove missing elements

        %-- Extra processing to account for partial elements -----%
        [dr,dr1,dr2] = obj.dr;
        switch dim % determine which dimension to sum over
            case 2 % integrate in column direction
                dr1 = obj.full2partial(dr1(:));
                dr0 = dr./dr1; % ~dr2
            case 1 % integrate in row direction
                dr2 = obj.full2partial(dr2(:));
                dr0 = dr./dr2; % ~dr1
        end
        C1 = bsxfun(@times,C1,dr0');
        
    end
    %=====================================================================%
    
    
    
    %== RESHAPE ======================================================%
    function x = reshape(obj, x)
    % RESHAPE  A simple function to reshape a vector based on the grid.
    %  NOTE: If the grid is partial, missing grid points are 
    %   filled with zeros. 
        
        x = obj.partial2full(x);
        
        x = reshape@Grid(obj, x);
        
    end
    %=================================================================%
    
    
    
    %== PLOT2D =======================================================%
    function [h, x] = plot2d(obj, x, f_tr, f_contf)
    % PLOT2D  Plots x as a 2D function on the grid.
    %  AUTHOR: Timothy Sipkens, 2018-11-21
        
        if ~exist('f_contf', 'var'); f_contf = []; end
        if ~exist('f_tr', 'var'); f_tr = []; end
        
        % Use default method to start plot.
        [h, x] = plot2d@Grid(obj, x, f_tr, f_contf);

        % Default for order of indices. 
        ptl = [0,obj.cut(1)];
        slope = obj.cut(2);
        if length(obj.cut)>2
            ptu = [0,obj.cut(3)];
            slope = [slope, obj.cut(4)];
        end
        
        % Transpose if flagged.
        if f_tr
            ptl = [obj.cut(1),0];
            if length(obj.cut)>2; ptu = [obj.cut(3),0]; end
            slope = 1 ./ slope;
        end
        
        % Add lines marking the edges of the partial grid.
        hold on;
        tools.overlay_line(obj, ptl, slope(1), ...
            'Color', [0.5,0.5,0.5]); % overlay partial grid limits, gray lines

         % If also a bottom cut.
        if length(obj.cut)>2
            tools.overlay_line(obj, ptu, slope(2), ...
                'Color', [0.5,0.5,0.5]); % add a gray line
        end
        hold off;
        
        if nargout>0; h = gca; end
        
    end
    %=================================================================%

    
    
    %== PARTIAL2FULL =================================================%
    function x_full = partial2full(obj ,x)
    % PARTIAL2FULL  Convert x defined on a partial grid to the full grid equivalent. 
    %  Uses zeros to fill the removed grid points.
        x_full = zeros(prod(obj.ne),1);
        t0 = setdiff((1:prod(obj.ne))', obj.missing);
        x_full(t0) = x;
    end
    %=================================================================%
    
    
    
    %== FULL2PARTIAL =================================================%
    function x = full2partial(obj, x)
    % FULL2PARTIAL  Convert x defined on a full grid to the partial grid equivalent. 
    %  Removes entries for missing indices.
        x(obj.missing,:) = [];
    end
    %=================================================================%
    
end

end

