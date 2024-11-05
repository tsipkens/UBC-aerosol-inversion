
% GEN_PMA_SP2_GRID  Evaluate kernel/transfer functions for PMA-SP2 to find A, exploiting grid structure.
% 
%  A = kernel.gen_pma_sp2_grid(GRID_B,GRID_I) uses the data grid, GRID_B,
%  and integration grid, GRID_I, which is generally higher resolution. 
% 
%  A = kernel.gen_pma_sp2_grid(GRID_B,GRID_I,PROP_PMA) uses the 
%  pre-computed particle mass analyzer properties specified by PROP_PMA. 
%  This structure can be generated manually or with the help of the 
%  kernel.prop_pma(...) function. 
% 
%  A = kernel.gen_pma_sp2_grid(...,'name',value) specifies a name-value
%  pairs that is passed on as input to the get_setpoint(...) function
%  from the tfer_pma folder. This specifies one quantity other than the
%  setpoint mass that is used to constrain the PMA operating point. For
%  example, {'Rm',10} specifies a resolution PMA resolution of 10, which is
%  used for all of the setpoints. 
%  If not provided, the function uses the default in the get_setpoint(...)
%  function in the tfer_pma folder.
% 
%  [A,SP] = kernel.gen_pma_sp2_grid(...) also outputs the PMA setpoint
%  structure (associated with the mat_tfer_pma  or tfer_pma folder) for 
%  the given grid.
%  
%  ------------------------------------------------------------------------
% 
%  NOTES:
%  1. This function exploits the grid structure to minimize the number of
%  transfer function evaluations. 
%  2. Cell arrays are used for Omega_mat and Lambda_mat in order to 
%  allow for the use of sparse matrices, which is necessary to 
%  store information on higher resolutions grids
%  (such as those used for phantoms).
%  3. By default, this function uses the analytical PMA transfer function 
%  corresponding to Case 1C from Sipkens et al. (Aerosol Sci. Technol. 2020b).
%  
%  AUTHOR: Timothy Sipkens, 2018-11-27

function [A,sp] = gen_pma_sp2_grid(grid_b, grid_i, prop_p, varargin)


%-- Parse inputs ---------------------------------------------------------%
% If not given, import default properties of PMA, 
% as selected by prop_pma function.
addpath tfer tfer/tfer-pma; % add mat-tfer-pma package to MATLAB path
if ~exist('prop_p','var'); prop_p = []; end
if isempty(prop_p); prop_p = prop_pma(); end
    
if or(isempty(varargin), length(varargin)~=2) % parse extra information for PMA
    error('Invalid additional information for PMA setpoint.');
end
%-------------------------------------------------------------------------%

    
%-- Parse measurement set points (b) -------------------------------------%
r_star = grid_b.elements;
m_star = r_star(:,2);
mrbc_star = r_star(:,1);
n_b = grid_b.ne;
N_b = prod(n_b); % length of data vector


%-- Generate grid for intergration ---------------------------------------%
n_i = grid_i.ne;
N_i = grid_i.Ne; % length of integration vector

r = grid_i.elements;
m = r(:,2);
mrbc = r(:,1);
d = (m.*1e-18./prop_p.m0).^...
    (1/prop_p.Dm);
    % invoke mass-mobility relation


%-- Start evaluate kernel ------------------------------------------------%
tools.textheader('Computing PMA-SP2 kernel');

%== Evaluate particle charging fractions =================================%
z_vec = (0:3)';
f_z = sparse(charger(d, z_vec)); % get fraction charged for d
n_z = length(z_vec);


%== STEP 1: Evaluate SP2 transfer function ===============================%
%   Note: The SP2 contribution is 1D and does not depend on the charge
%   state. It is boxcar function that takes into account discretization
%   only. 
disp(' Computing SP2 contribution...');
tools.textbar([0,n_b(1)]); % initiate textbar
Omega_mat = sparse(n_b(1),n_i(1));% pre-allocate for speed
for ii=1:n_b(1)
    Omega_mat(ii,:) = max(...
        min(grid_i.nodes{1}(2:end),grid_b.nodes{1}(ii+1))-... % lower bound
        max(grid_i.nodes{1}(1:(end-1)),grid_b.nodes{1}(ii))... % upper bound
        ,0)./(grid_b.nodes{1}(ii+1)-grid_b.nodes{1}(ii)); % normalize by SP2 bin size
    tools.textbar([ii,n_b(1)]);
end

[~,jj] = max(mrbc==grid_i.edges{1},[],2);
Omega_mat = Omega_mat(:,jj);
    % repeat transfer function for repeated mass in grid_i

disp(' Completed SP2 contribution.');
disp(' ');
%=========================================================================%


%== STEP 2: Evaluate PMA transfer function ===============================%
disp(' Computing PMA contribution:');

tools.textbar([0, n_z]); % initiate textbar
Lambda_mat = cell(1, n_z); % pre-allocate for speed, one cell entry per charge state
sp = get_setpoint(prop_p,...
    'm_star', grid_b.edges{2} .* 1e-18, varargin{:}); % get PMA setpoints

for kk=1:n_z % loop over the charge state
    Lambda_mat{kk} = tfer_pma(...
        sp, m, d, z_vec(kk), prop_p);
            % PMA transfer function

    tools.textbar([kk, n_z]);
end
disp(' Complete.');
disp(' ');

%{
% Diagnostic plot. 
tmp = Lambda_mat{1};
for ii=2:length(Lambda_mat)
    tmp = tmp + Lambda_mat{ii};
end
imagesc(tmp);
%}
%=========================================================================%


%== SETP 3: Combine to compile kernel ====================================%
disp(' Compiling kernel ...');
K = sparse(N_b,N_i);
for kk=1:n_z
    [~,i1] = max(m_star==grid_b.edges{2},[],2); % index corresponding to PMA setpoint
    [~,i2] = max(mrbc_star==grid_b.edges{1},[],2); % index correspondng to SP2 setpoint
    
    K = K+f_z(kk,:).*... % charging contribution
        Lambda_mat{kk}(i1,:).*... % PMA contribution
        Omega_mat(i2,:); % SP2 contribution
end
tools.textdone();  % print orange DONE
%=========================================================================%


dr_log = grid_i.dr; % area of integral elements in [logm,logd]T space
A = bsxfun(@times,K,dr_log'); % multiply kernel by element area
A = sparse(A); % exploit sparse structure

tools.textheader;


end



