
% BUILD  A general function to compute a 2D transfer function. 
%  
%  NOTE: This function generally requires that the type of particle size is
%  specified in grid.type = {}. The default if not specified is to assign
%  grid.type = {'dm', 'mp'}, corresponding to a mass-mobility scenario.
%  Other cases require identifying the dimensions of the grid. 
%  
%  A = kernel.build_grid(grid_i, z_vec, ...) builds a kernel based
%  on the grid for evaluation (GRID_I). Setpoints should be specified 
%  for each classifer Charge states to evaluate at are provided in Z_VEC. 
%  When Z_VEC is empty, the function defaults to Z_VEC = 1:3. Remaining 
%  elements are required and consistitute name-cell pairs for each classifier. 
%  
%  [A,AC] = kernel.build_grid(...) adds an output AC that is not summed
%  ovre the charge states, which is useful for advanced analysis. 
%  
%  A = kernel.build_grid(..., 'pma', {M_STAR, PROP_P, VARARGIN})
%  Builds a PMA contribution to the kernel using the properties in PROP_P
%  and a set of secondary information for the setpoints 
%  (e.g., {M_STAR, PROP_P, 'Rm', Rm}). M_STAR is expected in fg. 
%  
%  A = kernel.build_grid(..., 'dma', {D_STAR, PROP_D, VARARGIN})
%  Builds a DMA contribution to the kernel using the properties in PROP_D
%  and a set of secondary information for the classifer evaluations in
%  VARARGIN that is passed directly to tfer_dma. D_STAR is expected in nm.
%  
%  A = kernel.build_grid(..., 'charger', {VARARGIN})
%  Builds a charger contribution to the transfer function. An empty cell
%  will use the default call to charger(...). This addition is required for
%  all classifiers that require charging for classification. 
%  This does not require a mobility diameter setpoint. 
%  
%  Other classifiers follow this template. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2024-01-11

function [A, Ac] = build(grid_i, z_vec, varargin)

if mod(length(varargin), 2) ~= 0; error('Wrong number of inputs.'); end

if ~exist('z_vec', 'var'); z_vec = []; end
if isempty(z_vec); z_vec = 1:3; end  % default charge states for evaluation

tools.textheader('Computing kernel')

nc = (length(varargin)/2);  % number of classifiers
Lambda{nc} = [];  % initialize Lambda_ii

% Set default indices, if type not specified (mass-mobility grid).
if isempty(grid_i.type)
    grid_i.type = {'mp', 'dm'};
end

dm_idx = find(strcmp(grid_i.type, 'dm'));
mp_idx = find(strcmp(grid_i.type, 'mp'));
da_idx = find(strcmp(grid_i.type, 'da'));

% Get single particle mass (incl. from density if supplied).
% Convert effective density and mobility diameter to mass.
if any(strcmp(grid_i.type, 'rho'))
    rho_idx = find(strcmp(grid_i.type, 'rho'));
    if ~any(strcmp(grid_i.type, 'mp'))
        m = (pi/6) .* grid_i.elements(:, rho_idx) .* ...
            grid_i.elements(:, dm_idx) .^ 3 * 1e-9;
    end

% Convert dynamic shape factor and mobility diameter to mass. 
elseif any(strcmp(grid_i.type, 'chi'))
    idx_p = find(strcmp(varargin, 'pma')) + 1;  % first find PMA inputs
    prop_p = varargin{idx_p}{2};

    chi_idx = find(strcmp(grid_i.type, 'chi'));
    if ~any(strcmp(grid_i.type, 'mp'))
        m = 1e18 .* dm_chi2mp(grid_i.elements(:, dm_idx) .* 1e-9, ...
            grid_i.elements(:, chi_idx), prop_p.rhom);
    end

elseif any(strcmp(grid_i.type, 'mp'))
    m = grid_i.elements(:, mp_idx);
end

% Handle if mobility diameter is not given directly (req'd for charging/PMA).
% Compute using known relationships.
if isempty(dm_idx)
    addpath autils;  % add aerosol utilities (autils) package
    
    % OPTION 1: Use da and mp to directly compute dm.
    if and(~isempty(da_idx), ~isempty(mp_idx))
        m = grid_i.elements(:, mp_idx);  % get mass from relevant dimension
        da = grid_i.elements(:, da_idx);  % get da from relevant dimension
        dm = mp_da2dm(m, da);  % fully constrained calculation
        dm2 = dm';
    
    % OPTION 2: Apply assumption of a mass-mobility relationship, 
    % which will be less precise. Necessary for PMA-SP2.
    elseif ~isempty(mp_idx)
        idx_p = find(strcmp(varargin, 'pma')) + 1;  % first find PMA inputs
        prop_p = varargin{idx_p}{2};  % extract prop_pma from PMA input
        
        % Then convert using mass-mobility relationship. 
        disp(' Invoking mass-mobility relationship to determine dm.');
        dm = mp2dm(m .* 1e-18, prop_p) .* 1e9;
        dm2 = dm';
    end
end


% Loop over the various classifiers. 
for ii=1:nc
    jj = 2*ii - 1;  % index of classifier

    switch varargin{jj}
        
        %== CHARGER ======================================================%
        case 'charger'
            disp(' Computing charger contribution ...');
            
            % Assign inputs.
            if isempty(dm_idx)  % then PMA without DMA
                dm = dm';  % inherit from previous PMA calc.
                dm2 = dm';
            else
                dm = grid_i.edges{dm_idx};
                dm2 = grid_i.elements(:, dm_idx);
            end

            f_z = charger(dm', z_vec, varargin{jj+1}{1:end}); % get fraction charged for d vector
            Lambda{ii} = permute(f_z, [3, 2, 1]);

            % Duplicate over other grid dimensions.
            [~,kk] = max(dm == dm2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            tools.textdone();


        %== SMPS =========================================================%
        case {'dma', 'smps'}
            disp(' Computing DMA contribution ...');
            
            % Unpack inputs.
            d_star2 = varargin{jj+1}{1};
            d_star = unique(d_star2)';
            d = grid_i.edges{dm_idx};  % points for integration
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_dma(d_star, d', z_vec, varargin{jj+1}{2:end});

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, dm_idx);
            [~,kk] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            [~,kk] = max(d_star == d_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();

        
        %== PMA ==========================================================%
        %   Currently assumes other dimension is mobility diameter. 
        case 'pma'
            disp(' Computing PMA contribution ...');

            % Unpack inputs.
            m_star = varargin{jj+1}{1};  % don't use unique(), as resolution may change
            prop_p = varargin{jj+1}{2};  % DMA properties
            
            % Handle mobility diameter.
            if ~isempty(dm_idx)  % use corresponding dimension of grid
                dm = grid_i.elements(:, dm_idx);
            else  % then likely PMA without DMA
                dm = (m .* 1e-18 ./ prop_p.rho0) .^ ...
                    (1/prop_p.Dm) .* 1e9;  % use mass-mobility
                disp('  Invoking mass-mobility relation for PMA.')
            end

            addpath 'tfer\tfer-pma';  % added to calculate sp
            sp = get_setpoint(prop_p,...  % get PMA setpoints
                'm_star', m_star .* 1e-18, ...  % mass from the grid
                varargin{jj+1}{3:end});  % extra name-value pair to specify setpoint
            
            % Find unique setpoints. 
            [spu, ~, kk] = unique([[sp.m_star]', [sp.Rm]'], 'rows');
            spu = get_setpoint(prop_p, 'm_star', spu(:,1), 'Rm', spu(:,2));
            
            Lambda{ii} = tfer_pma(...
                spu, m, dm, z_vec, prop_p);

            % Duplicate over repeat entries.
            Lambda{ii} = Lambda{ii}(kk,:,:);

            tools.textdone();
            

        %== AAC ==========================================================%
        case 'aac'
            disp(' Computing AAC contribution ...');
            
            addpath 'autils';
            
            % Assign inputs.
            d_star2 = varargin{jj+1}{1};

            % Copy over prop structure.
            prop = varargin{jj+1}{2};
            f = {'Qa', 'Qs', 'Qsh', 'Qexh'};
            nf = length(f);
            
            % Get unique flow/da_star combinations.
            % Build unique vector.
            u = d_star2;
            for ff=1:nf  % extend flows to allow for vector input
                u = [u, prop.(f{ff}) .* ones(size(d_star2))];
            end
            [u, ~, kk] = unique(u, 'rows');  % get unique rows
            d_star = u(:, 1)';
            for ff=1:nf  % deconstruct flows
                prop.(f{ff}) = u(:, 1 + ff)';
            end
            
            d = grid_i.edges{da_idx};  % points for integration
            
            % Evaluate transfer function.
            Lambda{ii} = tfer_aac(d_star, d', prop, varargin{jj+1}{3:end});
            Lambda{ii} = permute(Lambda{ii}, [2,1]);

            % Duplicate over other grid dimensions.
            d2 = grid_i.elements(:, da_idx);
            [~,ll] = max(d == d2, [], 2);
            Lambda{ii} = Lambda{ii}(:,ll,:);
            
            % [~,kk] = max(d_star == d_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);
            
            tools.textdone();
            
        
        %== BIN ==========================================================%
        %   When data input is binned (e.g., SP2 data).
        case {'bin', 'sp2'}
            disp(' Computing binned contribution ...');

            % Unpack inputs.
            s_star2 = varargin{jj+1}{1};
            s_star = unique(s_star2)';

            s_idx = varargin{jj+1}{2};
            s = grid_i.edges{s_idx};  % points for integration
            
            Lambda{ii} = full(tfer_bin(s_star', s'));
            
            % Duplicate over other grid dimensions.
            s2 = grid_i.elements(:, s_idx);
            [~,kk] = max(s == s2, [], 2);
            Lambda{ii} = Lambda{ii}(:,kk,:);
            
            % Copy to identical data points.
            [~,kk] = max(s_star == s_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);
            
            % Avoid double counting bin width during later ".* dr".
            [~, dr1, dr2] = grid_i.dr;
            if s_idx == 1; dr = dr1(:);
            else dr = dr2(:); end

            if isa(grid_i, 'PartialGrid')
                dr = grid_i.full2partial(dr);
            end
            Lambda{ii} = Lambda{ii} ./ dr';

            tools.textdone();


        case {'sp2-frBC'}
            disp(' Computing SP2 (frBC) contribution ...');

            % Unpack inputs.
            s_star2 = varargin{jj+1}{1};
            s_star = unique(s_star2)';

            s_idx = find(strcmp(grid_i.type, 'frBC'));
            alt_idx = 3 - s_idx;

            frBC = grid_i.edges{s_idx};
            Lambda{ii} = zeros(length(s_star), grid_i.Ne);
            for ss=1:length(grid_i.edges{alt_idx})
                mp = grid_i.edges{alt_idx}(ss);
                mrBC = frBC .* mp;
                
                Lam_ss = full(tfer_bin(s_star', mrBC'));

                frBC2 = grid_i.elements(:, s_idx);
                mp2 = grid_i.elements(:, alt_idx);
                [~,kk] = max(and(frBC == frBC2, mp == mp2), [], 2);
                Lam_ss = Lam_ss(:,kk);

                % Remove null rows.
                Lam_ss(:, ~any(and(frBC == frBC2, mp == mp2), 2)) = 0;

                Lambda{ii} = Lambda{ii} + Lam_ss;
            end
            
            % Copy to identical data points.
            [~,kk] = max(s_star == s_star2, [], 2);
            Lambda{ii} = Lambda{ii}(kk,:,:);
            
            % Avoid double counting bin width during later ".* dr".
            [~, dr1, dr2] = grid_i.dr;
            if s_idx == 1; dr = dr1(:);
            else dr = dr2(:); end

            if isa(grid_i, 'PartialGrid')
                dr = grid_i.full2partial(dr);
            end
            Lambda{ii} = Lambda{ii} ./ dr';

            tools.textdone();

    end
end

% Loop over the various classifiers again to compile kernel.
disp(' Compiling kernel ...')
A = Lambda{1};  % initialize with first contribution
for ii=2:nc  % loop over other contributions
    A = A .* Lambda{ii};
end

% If second output selected, provide output prior to summing over charge.
if nargout > 1
    Ac = A;
end

A = sum(A, 3);  % sum over charge states

A = A .* grid_i.dr';  % multiply kernel by element area
A = sparse(A);  % exploit sparse structure in subsequent calculations

tools.textdone();
tools.textheader();

end



function [mp, dve] = dm_chi2mp(dm, chi, rhom)
% DM_CHI2MP  Use mobility diameter and shape factor to get volume-eq. diameter.
%  AUTHOR: Timothy Sipkens, 2024-10-24

addpath autils;

dve = dm ./ chi;

for ii=1:length(dve)
    fun = @(dve) dve - dm(ii) / chi(ii) .* Cc(dve) ./ Cc(dm(ii));
    dve(ii) = fzero(fun, dve(ii));
end

mp = rhom .* (pi/6) .* dve .^ 3;

end


