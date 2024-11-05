
% PHANTOM  Creates a bivariate phantom (typically bivariate lognormal).
%  Can be used to represent multimodal phantoms.
% 
%  Use the Phantom.eval(grid) method to evalute the phantom on a 
%  specific Grid.
%  
%  PHA = Phantom(TYPE) generates a phantom distribution, where type is one
%  of multiple options below, specifying different ways of building the
%  phantom. 
%  
%  PHA = Phantom(TYPE, GRID, ...) adds a instance of the Grid class on which
%  that phantom will be evaluated. If GRID is not given, the phantom is
%  simply not evaluated at the time of contrusction (but can be evaluated
%  on a specified grid using Phantom.eval(grid). 
%  
%  PHA = Phantom(TYPE, GRID, ~, ~, W) constructs a phantom with mode 
%  weights specified by W, which is an array with one weight per mode. 
%  By default, all modes are given equal weight and selected such that the 
%  sum of the weights is unity. 
%  
%  ------------------------------------------------------------------------
%  
%  OPTION 1: The 'standard' parameterization
% 
%  PHA = Phantom('standard', [], MU, SIGMA) creates a phantom with a mean of 
%  MU and covariance of SIGMA. This is the most general definition. In the 
%  standard representation, all modes are taken as bivariate lognormal, 
%  where MU is log10(...) of the geometric means and SIGMA is the covariance 
%  in logspace. Bimodal distributions are given by stacking the means and
%  covariance in the next dimension. For MU, this results in one row per
%  set of means. For SIGMA, this involves concatenating covariance matrices
%  in the third dimension, e.g., SIGMA = cat(3, SIGMA1, SIGMA2). Examples: 
%  
%  ```
%  pha = Phantom('standard', [], [-0.105,1.70], [0.205,0.0641; 0.0641,0.0214])
%  pha = Phantom('standard', [], ...
%  	[0,2; 1,2.5],... % distribution means
%  	cat(3, [0.205,0.0641; 0.0641,0.0214],... % covaraince of mode no. 1
%  	[0.126,0.0491; 0.0491,0.0214])); % covaraince of mode no. 2
%  ```
%  
%  See Sipkens et al., J. Aerosol Sci. (2020) for more information
%  on the bivariate lognormal distribution. 
%  
%  ------------------------------------------------------------------------
%  
%  OPTION 2: The 'mass-mobility' parameterization
%  
%  PHA = Phantom('mass-mobility', [], P, MODES) creates a mass-mobility 
%  phantom using a P data structure. MODES is a cell array, with entries
%  specifying the type of distribution the P struct represents for  
%  each mode. Specifically: 
%   - 'logn' for a bivariate lognormal mode or 
%   - 'norm' fpr conditionally-normal distribution, where the mobility 
%     diameter distribution is lognormal and the conditional mass  
%     distribution is normal (as per Buckley et al., 2017). 
%  
%  This is built specifically for mass-mobility distributions. 
%  The required fields for the P structure are: 
%  
%  1. `dg` =  Mean mobility diameter  
%  2. `sg` =  Standard deviation of the mobility diameter
%  3. `Dm` =  Mass-mobility exponent
%  4. *Either:*
%   `sm`   =  Standard deviation of the particle mass
%   `smd`  =  Standard deviation of the conditional mass distribution  
%  5.  *Either:*
%   `mg`   =  Mean particle mass
%   `rhog` =  Effective density of at the mean mobility diameter
%  
%  For lognormal modes, means should be geometric means and standard 
%  deviations should be geometric standard deviations. Remaining entries 
%  of the `p` structure will be filled using the `Phantom.fill_p(...)` 
%  method. (We note that `p = Phantom.fill_p(p);` can be used to fill out 
%  the `p` structure without the need to create an instance of the Phantom 
%  class.) Examples:
%  
%  ```
%  p.dg = 50; p.rhog = 12000; % geometric means
%  p.sg = 1.4; p.smd = 1.3; % geometric standard deviations
%  p.Dm = 3; % mass-mobility exponent
%  pha = Phantom('mass-mobility',[],p,{'logn'}); % create phantom
%  
%  % For Phantom '1'":
%  % mode 1           mode 2
%  p(1).dg = 50;      p(2).dg = 200;
%  p(1).rhog = 12000; p(2).rhog = 500;
%  p(1).sg = 1.4;     p(2).sg = 1.4;
%  p(1).smd = 1.3;    p(2).smd = 1.3;
%  p(1).Dm = 3;       p(2).Dm = 2.3;
%  pha = Phantom('mass-mobility',[],p,{'logn','logn'});
%  ```
%  
%  ------------------------------------------------------------------------
%  
%  OPTION 3: Preset phantoms
% 
%  PHA = Phantom(NAME) generates a preset phantom object specified by 
%  NAME, a string that is typically an integer. For example, '1' generates 
%  Phantom 1 from Sipkens et al., J. Aerosol Sci. (2020). One can also
%  supply integers in the range of [1,4] for the phantoms in the
%  aforementioned paper. The corresponding function are loaded using an
%  external `presets` function. 
%  
%  ```
%  phantom = Phantom('1');
%  ```
%  
%  ------------------------------------------------------------------------
%  
%  OPTION 4: Using the Phantom class's fit methods
%  
%  For experimental data, the Phantom class can also be used to derive 
%  morphological parameters from the reconstructions. 
%  
%  Of particular note, the `Phantom.fit(...)` method, which is defined 
%  external to the main definition of the Phantom class, takes a 
%  reconstruction, `x` and the grid on which it is defined and creates a 
%  bivariate lognormal phantom that most resembles the data. This done 
%  using least squares analysis. The `p` structure of the Phantom class 
%  then contains many of the morphological parameters of interest to 
%  practitioners measuring mass-mobility distributions. 
%  
%  The `Phantom.fit2(...)` method can be used in an attempt to derive 
%  multimodal phantoms for the data. This task is often challenging, such 
%  that the method may need tuning in order to get distributions that 
%  appropriately resemble the data. 
%  
%  ------------------------------------------------------------------------
%  
%  Converting between the 'standard' and 'mass-mobility' parameterizations
%  
%  In both cases, creating an instance of the class will also contain the 
%  information corresponding to the other creation method (e.g. using a 
%  `p` structure, the class constructor will determined the corresponding 
%  mean and covariance information and store this in `Phantom.mu` and 
%  `Phantom.Sigma`). This can be demonstrated by investigating the examples 
%  provided in the proceeding sections, which generate the same phantom, 
%  save for rounding errors.  Conversion between the '**standard**' 
%  parameterization and the '**mass-mobility**' parameterizations can be 
%  accomplished using the `Phantom.cov2p(...)` method of the Phantom class 
%  and vice versa using the `Phantom.p2cov(...)` method of the Phantom 
%  class. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2019-07-08


classdef Phantom

%-- Phantom properties -----------------------------------------------%
properties
    type = [];      % optional name for the phantom
    modes = {};     % types of distribution for each mode, e.g. {'logn','logn'}
    n_modes = [];   % number of modes
    w = [];         % weighting for each mode

    mu = [];        % center of bivariate distribution
    Sigma = [];     % covariance of bivariate distribution
    R = [];         % correlation matrix

    x = [];         % phantom evaluated on default grid
    grid = [];      % default grid the phantom is to be represented
                    % on generally a high resolution instance of the
                    % Grid class
    
    % Parameters relevant to mass-mobility distributions
    % fields include mg, sm, sd, rho_100, etc.
    p = struct();
    
end



methods
    %== PHANTOM ======================================================%
    function [obj] = Phantom(type_name, span_grid, mu_p, Sigma_modes, w)
        
        %-- Parse inputs ---------------------------------------------%
        if nargin==0; return; end  % return empty phantom
        
        if ~exist('span_grid','var'); span_grid = []; end
        
        if ~exist('w','var'); w = []; end
        
        % If type is a number, cover to string for future comparisons.
        if isnumeric(type_name); type_name = num2str(type_name); end
        %-------------------------------------------------------------%
        
        %== Assign parameter values - 3 options ======================%
        switch type_name
            
            %-- OPTION 1: Standard bivariate lognormal distribution --%
            case {'standard'}
                n_modes = size(mu_p,1);
                
                obj.mu = mu_p;
                obj.Sigma = Sigma_modes;
                obj.R = obj.cov2corr(obj.Sigma);
                
                for ii=1:n_modes; obj.modes{ii} = 'logn'; end
                
                p = obj.cov2p(obj.mu,obj.Sigma,obj.modes);
                    % mass-mobility equivlanet params.
                obj.p = obj.fill_p(p); % get mg as well
                
                
            %-- OPTION 2: Using a mass-mobility parameter set (p) ----%
            case {'mass-mobility'} % for custom mass-mobility phantom
                                   % specified using a p structure
                obj.type = 'mass-mobility';
                obj.modes = Sigma_modes;
                obj.p = obj.fill_p(mu_p); % fill out p structure
                
                if ~any(strcmp('cond-norm',Sigma_modes))
                    [obj.mu,obj.Sigma] = obj.p2cov(obj.p,obj.modes);
                    obj.R = obj.cov2corr(obj.Sigma);
                end
                
                
            %-- OPTION 3: Use a preset or sample distribution --------%
            otherwise % check if type is a preset phantom
                [p, modes, type_name, w] = obj.presets(type_name);
                if isempty(p); error('Invalid phantom call.'); end
                
                obj.type = type_name;
                obj.modes = modes;
                
                obj.p = obj.fill_p(p);
                
                if ~any(strcmp('cond-norm',modes))
                    [obj.mu,obj.Sigma] = obj.p2cov(obj.p,obj.modes);
                    obj.R = obj.cov2corr(obj.Sigma);
                end
        end
        
        obj.n_modes = length(obj.modes); % get number of modes
        
        if isempty(w) % assign mode weightings
            obj.w = ones(obj.n_modes, 1) ./ obj.n_modes; % evely distribute modes
        else
            obj.w = w./sum(w); % normalize weights and assign
        end
            
        
        %-- Generate a grid to evaluate phantom on -------------------%
        if isa(span_grid, 'Grid') % if grid is specified
            obj.grid = span_grid;
        elseif ~isempty(span_grid) % if span is specified, create grid
            n_t = [540,550]; % resolution of phantom distribution
            obj.grid = Grid(span_grid, ...
                n_t, 'logarithmic'); % generate grid of which to represent phantom
        end


        %-- Evaluate phantom -----------------------------------------%
        if ~isempty(span_grid)
            if any(strcmp('cond-norm', obj.modes))
                obj.x = obj.eval_p(obj.p);
                    % special evaluation for conditional normal conditions
            else
                obj.x = obj.eval;
            end
        end
    end
    %=================================================================%



    %== MG_FUN =======================================================%
    %   Function to evaluate mg as a function of mobility diameter.
    %   Calculation is based on the empirical mass-mobility relation.
    %   Author:  Timothy Sipkens, 2019-07-19
    function mg = mg_fun(obj,d)

        d_size = size(d);
        if d_size(2)>d_size(1); d = d'; d_size = size(d); end
            % transpose the mobility if necessary

        mg = zeros(d_size(1),obj.n_modes);
        rho = obj.rho_fun(d);
        for ll=1:obj.n_modes
            mg(:,ll) = 1e-9.*rho(:,ll).*pi./6.*(d.^3);
                % output in fg
        end

    end
    %=================================================================%



    %== RHO_FUN ======================================================%
    %   Function to evaluate the effective density as a function of
    %   mobility diameter.
    %   Author:  Timothy Sipkens, 2019-07-19
    function rho = rho_fun(obj,d)

        rho = zeros(length(d),obj.n_modes);
        for ll=1:obj.n_modes
            rho(:,ll) = 6*obj.p(ll).k./(pi.*d.^(3-obj.p(ll).Dm));
                % output in kg/m3
        end

    end
    %=================================================================%



    %== PLOT =========================================================%
    %   Plots the phantom distribution.
    %   Author:     Timothy Sipkens, 2019-07-08
    function [h] = plot(obj)
        h = obj.grid.plot2d_marg(obj.x);
    end
    %=================================================================%



    %== EVAL =========================================================%
    %   Generates a distribution from the phantom mean and covariance.
    %   Author:  Timothy Sipkens, 2019-10-29
    %
    %   Note: This method does not work for conditionally-normal distributions
    %       (which cannot be defined with mu and Sigma).
    function [x] = eval(obj, grid_vec, w)
        
        %-- Parse inputs ---------------------------------------------%
        if ~exist('w','var'); w = []; end
        if isempty(w); w = obj.w; end
        if isempty(w); w = ones(obj.n_modes,1)./obj.n_modes; end
            % weight modes evenly
            
        if ~exist('grid_vec','var'); grid_vec = []; end
        if isempty(grid_vec); grid_vec = obj.grid; end % use phantom grid
        if isempty(grid_vec); error('For an empty Phantom, grid_vec is required.'); end
            % if an empty phantom (e.g. Phantom.eval_p(p);)
        
        if isa(grid_vec,'Grid'); vec = grid_vec.elements; % get element centers from grid
        else; vec = grid_vec; % if a set of element centers was provided directly
        end
        
        mu0 = obj.mu;
        Sigma0 = obj.Sigma;
        %-------------------------------------------------------------%
        
        if ~iscell(obj.mu); mu0 = {obj.mu};
        else; mu0 = obj.mu;
        end
        
        m_vec = vec(:,1); % element centers in mass
        d_vec = vec(:,2); % element centers in mobility

        %-- Assign other parameters of distribution ------------------%
        x = zeros(size(m_vec));
        for ll=1:obj.n_modes % loop through distribution modes
            x = x + w(ll).*... % add new mode, reweighting accordingly
                mvnpdf(log10([m_vec,d_vec]),...
                obj.mu(ll,:),...
                obj.Sigma(:,:,ll));
        end
    end
    %=================================================================%



    %== EVAL_P =======================================================%
    %   Generates a distribution from p as required for conditionally-
    %   normal modes.
    %   Author:  Timothy Sipkens, 2019-10-29
    function [x] = eval_p(obj,p,grid_vec,w)
        
        %-- Parse inputs ---------------------------------------------%
        if ~exist('w','var'); w = []; end
        if isempty(w); w = obj.w; end
        if isempty(w); w = ones(obj.n_modes,1)./obj.n_modes; end
            % weight modes evenly
        
        if ~exist('p','var'); p = []; end
        if isempty(p); p = obj.p; end % use p values from given phantom
        if isempty(p); error('For an empty Phantom, p is required.'); end 
            % if an empty phantom (e.g. Phantom.eval_p;)
        
        if ~exist('grid','var'); grid_vec = []; end
        if isempty(grid_vec); grid_vec = obj.grid; end % use phantom grid
        if isempty(grid_vec); error('For an empty Phantom, grid_vec is required.'); end
            % if an empty phantom (e.g. Phantom.eval_p(p);)
        
        if isa(grid_vec,'Grid'); vec = grid_vec.elements; % get element centers from grid
        else; vec = grid_vec; % if a set of element centers was provided directly
        end
        %-------------------------------------------------------------%
        
        m_vec = vec(:,1); % element centers in mass
        d_vec = vec(:,2); % element centers in mobility
        
        m_fun = @(d,ll) log(p(ll).m_100.*((d./100).^p(ll).Dm));
            % geometric mean mass in fg as a function of d
        
        %-- Evaluate phantom mass-mobility distribution ---------------%
        x = zeros(length(m_vec),1); % initialize distribution parameter
        for ll=1:obj.n_modes % loop over different modes
            if strcmp(obj.modes{ll},'logn')
                p_m = lognpdf(m_vec,m_fun(d_vec,ll),log(p(ll).smd));
            else
                p_m = normpdf(m_vec,...
                    exp(m_fun(d_vec,ll)),p(ll).smd.*...
                    exp(m_fun(d_vec,ll)));
            end
            
            x = x + ...
                w(ll).*p_m.*...
                lognpdf(d_vec,log(p(ll).dg),log(p(ll).sg));
        end
        
        %-- Reweight modes and transform to log-log space ------------%
        x = x.*(d_vec.*m_vec).*log(10).^2;
            % convert to [log10(m),log10(d)]T space
    end
    %=================================================================%
    
    
    
    %== PLUS =========================================================%
    %   Adds two phantoms.
    %   Author:  Timothy Sipkens, 2020-03-24
    function objn = plus(obj1,obj2,w)
        
        if ~exist('w','var'); w = []; end
        if isempty(w); w = [1,1]./2; else; w = w./sum(w); end
        w = [w(1).*obj1.w(:);w(2).*obj2.w(:)];
        
        span = [min(obj1.grid.span(:,1),obj2.grid.span(:,1)),...
                max(obj1.grid.span(:,2),obj2.grid.span(:,2))];
                % update span to incorporate both Phantoms
        
                
        if and(~isempty(obj1.mu),~isempty(obj2.mu)) % bivariate lognormal phantoms
            
            objn = Phantom('standard',span,...
                [obj1.mu;obj2.mu],...
                cat(3,obj1.Sigma,obj2.Sigma),...
                w);
            
            
        else % at least one mode is not bivariate lognormal
            objn = Phantom('mass-mobility',span,...
                [obj1.p(:);obj2.p(:)]',...
                [obj1.modes,obj2.modes],...
                w);
        end
    end
    
    
    
    %== MASS2RHO =====================================================%
    %   Convert a mass-mobility phantom to an effective density-mobility
    %   phanatom. Output is a new phantom in the transformed space.
    %   Author:  Timothy Sipkens, 2019-10-31
    function [phantom] = mass2rho(obj,grid_rho)

        A = [1,-3;0,1]; % corresponds to mass-mobility relation
        
        for ll=1:obj.n_modes
            mu_rhod(ll,:) = (A*obj.mu(ll,:)'+[log10(6/pi)+9;0])';
            Sigma_rhod(:,:,ll) = A*obj.Sigma(:,:,ll)*A';
        end
        
        phantom = Phantom('standard',grid_rho,mu_rhod,Sigma_rhod);

    end
    %=================================================================%
    
    
    
    %== RHO2MASS =====================================================%
    %   Convert an effective density-mobility phantom to a mass-mobility
    %   phanatom. Output is a new phantom in the transformed space.
    %   Author:  Timothy Sipkens, 2020-06-01
    function [phantom] = rho2mass(obj,grid_mm)

        A = [1,3;0,1]; % corresponds to mass-mobility relation
        
        for ll=1:obj.n_modes
            mu_mm(ll,:) = (A*obj.mu(ll,:)'-[log10(6/pi)+9;0])';
            Sigma_mm(:,:,ll) = A*obj.Sigma(:,:,ll)*A';
        end
        
        phantom = Phantom('standard',grid_mm,mu_mm,Sigma_mm);

    end
    %=================================================================%
end



methods (Static)
    %== PRESETS ======================================================%
    %   Returns a set of parameters for preset/sample phantoms.
    %   AUTHOR: Timothy Sipkens, 2019-10-31, 2024-07-12
    function [p, type, name, w] = presets(name)
    
    fname = ['@Phantom/json/', name, '.json'];
    
    if ~isfile(fname)
        error('The specified phantom config file does not exist.');
    end
    
    % Read in phantom information.
    s = io.read_json(fname);
    
    % Unpack phantom information. 
    if isfield(s, 'p')
        name = s.name;
        p = s.p;
        type = s.type;
    
        if isfield(s, 'w'); w = s.w; end
    
    elseif isfield(s, 'Gam')
        % Not currently available.
    end
    
    % If no weight, then assign as empty for output. 
    if ~exist('w', 'var'); w = []; end
    
    end
    %=================================================================%
    
    
    %== FIT (External definition) ====================================%
    %   Fits a phantom to a given set of data, x, defined on a given grid, 
    %       or vector of elements. Outputs a fit phantom object.
    [phantom,N,y_out,J] = fit(x,vec_grid,logr0,Sigma);
    %=================================================================%
    
    %== FIT2 (External definition) ===================================%
    %   Fits a multimodal phantom object to a given set of data, x, 
    %       defined on a given grid or vector of elements.
    %   Outputs a fit phantom object.
    [phantom,N,y_out,J] = fit2(x,vec_grid,n_modes,logr0);
    %=================================================================%
    
    %== FIT2_RHO (External definition) ===============================%
    %   Fits a multimodal phantom object to a given set of data, x, 
    %       defined on a given grid or vector of elements. 
    %   Outputs a fit phantom object.
    %   Tuned specifically for effective density-mobility distributions
    %       (e.g. for negative or nearly zero correlation)
    [phantom,N,y_out,J] = fit2_rho(x,vec_grid,n_modes,logr0);
    %=================================================================%
    
    %== FIT_GMM (External definition) ================================%
    %	Fits a phantom to a given set of data, x, defined on a given grid.
    %   Uses sampling and k-means to fit a Guassian mixture model.
    %   Outputs a fit phantom object.
    [phantom,N,s] = fit_gmm(x,grid,k);
    %=================================================================%
    
    

    %== FILL_P =======================================================%
    %   Generates the remainder of the components of p.
    %   Requires a minimum of: 
    %   Author:  Timothy Sipkens, 2019-10-30
    function p = fill_p(p)
        
        n_modes = length(p);
        
        %-- Assign other parameters of distribution ------------------%
        for ll=1:n_modes % loop through distribution modes
            
            %-- Handle mg/rhog ------------%
            if ~isfield(p(ll),'rhog'); p(ll).rhog = []; end
            if isempty(p(ll).rhog)
                p(ll).rhog = p(ll).mg/(1e-9*pi/6*p(ll).dg^3);
            end
            p(ll).mg = 1e-9*p(ll).rhog*pi/6*...
                (p(ll).dg^3); % geometric mean mass in fg
            
            %-- Handle sm/smd -------------%
            if ~isfield(p(ll),'smd'); p(ll).smd = []; end
            if isempty(p(ll).smd)
                p(ll).smd = 10^sqrt(log10(p(ll).sm)^2 - ...
                    p(ll).Dm^2*log10(p(ll).sg)^2);
            end
            p(ll).sm = 10^sqrt(log10(p(ll).smd)^2+...
                p(ll).Dm^2*log10(p(ll).sg)^2);
            
            %-- Other parameters ----------%
            p(ll).rho_100 = p(ll).rhog*(100/p(ll).dg)^(p(ll).Dm-3);
            p(ll).m_100 = 1e-9*p(ll).rho_100*pi/6*100^3;
            p(ll).rhog = p(ll).rho_100*((p(ll).dg/100)^(p(ll).Dm-3));
            p(ll).k = p(ll).m_100/(100^p(ll).Dm);
            
            % Alternative variable names.
            p(ll).zet = p(ll).Dm;
            p(ll).m0 = p(ll).k;
        end
        
    end
    %=================================================================%



    %== VEC2P ========================================================%
    %   Function to format phantom parameters from a vector, t.
    %   Author:  Timothy Sipkens, 2019-07-18
    function [p] = vec2p(vec,n_modes)
        t_length = length(vec);
        n = n_modes*t_length;

        p.dg = vec(1:t_length:n);
        p.sg = vec(2:t_length:n);
        p.rho_100 = vec(3:t_length:n);
        p.smd = vec(4:t_length:n);
        p.Dm = vec(5:t_length:n);
    end
    %=================================================================%



    %== COV2P ========================================================%
    %   Function to convert covariance matrix and mean to p.
    %   Author:  Timothy Sipkens, 2019-10-29
    function [p,Dm,l1,l2] = cov2p(mu,Sigma,modes)

        if ~exist('modes','var'); modes = [];end
        if isempty(modes); modes = repmat({'logn'},[1,length(mu)]); end

        p = [];
        for ll=length(modes):-1:1 % loop through modes
            p(ll).dg = 10.^mu(ll,2);
            p(ll).mg = mu(ll,1);

            if strcmp(modes{ll},'logn') % if lognormal distribution, convert to geometric mean
                p(ll).mg = 10.^p(ll).mg;
            end

            p(ll).sg = 10^sqrt(Sigma(2,2,ll));
            p(ll).sm = 10^sqrt(Sigma(1,1,ll));

            R12 = Sigma(1,2,ll)/...
                sqrt(Sigma(1,1,ll)*Sigma(2,2,ll));
            p(ll).smd = 10^sqrt(Sigma(1,1,ll)*(1-R12^2));
                % conditional distribution width

            p(ll).Dm = Sigma(1,2,ll)/Sigma(2,2,ll);
                % corresponds to slope of "locus of vertical"
                % (Friendly, Monette, and Fox, 2013)
                % can be calculated as Dm = corr*sy/sx

            % t0 = eigs(rot90(Sigma(:,:,ll),2),1);
            % p(ll).ma = (t0-Sigma(2,2,ll))./Sigma(1,2,ll);
                % calculate the major axis slope

            p(ll).rhog = p(ll).mg/(pi*p(ll).dg^3/6)*1e9;
        end
        
        p = Phantom.fill_p(p);
        Dm = [p.Dm];
        l1 = log10([p.sm]);
        l2 = log10([p.sg]);
    end
    %=================================================================%



    %== P2COV ========================================================%
    %   Function to convert p to a covariance matrix and mean.
    %   Author:  Timothy Sipkens, 2019-10-29
    function [mu,Sigma] = p2cov(p,modes)

        for ll=length(p):-1:1
            mu(ll,:) = log10([p(ll).mg,p(ll).dg]);
                % use geometric mean

            if strcmp(modes{ll},'logn')
                % R12 = (1+1/(p(ll).Dm^2)*...
                %     (log10(p(ll).smd)/log10(p(ll).sg))...
                %     ^2)^(-1/2);
                Sigma(:,:,ll) = inv([(1/log10(p(ll).smd))^2,...
                    -p(ll).Dm/log10(p(ll).smd)^2;...
                    -p(ll).Dm/log10(p(ll).smd)^2,...
                    1/log10(p(ll).sg)^2+p(ll).Dm^2/log10(p(ll).smd)^2]);
            
            else  % for non-bivariate lognormal distributions, approximate
                Sigma(:,:,ll) = inv([(1/p(ll).smd)^2,...
                    -p(ll).Dm/p(ll).smd^2;...
                    -p(ll).Dm/p(ll).smd^2,...
                    1/log10(p(ll).sg)^2+p(ll).Dm^2/p(ll).smd^2]);
            end
        end
    end
    %=================================================================%



    %== COV2CORR =====================================================%
    %   Function to convert covariance matrix to correlation matrix.
    %   Author:  Timothy Sipkens, 2019-10-29
    function R = cov2corr(Sigma)
        for ll=length(Sigma(1,1,:)):-1:1
            R12 = Sigma(1,2,ll)/...
                sqrt(Sigma(1,1,ll)*Sigma(2,2,ll));
                % off-diagonal correlation

            R(:,:,ll) = diag([1,1])+rot90(diag([R12,R12]));
                % form correlation matrix
        end
    end
    %=================================================================%
    
end

end
