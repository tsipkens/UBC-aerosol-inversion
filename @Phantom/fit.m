
% FIT  Fits a phantom to a given set of data, x, defined on a given grid. 
%      Outputs a fit phantom object.
% 
% Inputs:
%   x           Input data (2D distribution data)
%   vec_grid    A `Grid` object on which input data is defined or a vector
%               of the coordiantes of the elements
%   logr0       Initial guess for center of each mode, expressed as a vector:
%               logr0 = [dim1_mode1,dim2_mode1]
%   Sigma       Estimate of the phantom covariance
%
% Outputs: 
%   pha         Output phantom object representing a fit bivariate lognormal
%   N           Scaling parameter that acts to scale the phantom to the data
%   y_out       Direct output from the fitting procedure
%   J           Jacobian of fitting procedure
%=============================================================%

function [pha, N, y_out, ci] = fit(x, vec_grid, logr0, Sigma)

tools.textheader('Fitting phantom object');

%-- Parse inputs ---------------------------------------------%
if isa(vec_grid,'Grid')
    grid = vec_grid;
    [~,vec1,vec2] = grid.vectorize();
else
    vec1 = vec_grid(:,1);
    vec2 = vec_grid(:,2);
    grid = [min(vec1),max(vec1);min(vec2),max(vec2)];
        % specify a span for a grid
end
if ~exist('logr0','var'); logr0 = []; end
if ~exist('Sigma','var'); Sigma = []; end
%-------------------------------------------------------------%


corr2cov = @(sigma,R) diag(sigma) * R * diag(sigma);

fun_pha = @(y) y(1) .* mvnpdf(...
    log10([vec1,vec2]), [y(2),y(3)],...
    corr2cov([y(4),y(5)], [1,sin(y(6));sin(y(6)),1]));
y0 = [max(x), ...  % C, scaling factor
    0, 2.3, ...  % log10(mg),log10(dg)
    0.6, 0.2, ...  % log10(sm),log10(sg)
    asin(0.99)];  % asin(corr), sin() helps to better span space

% Update centers of distributions, if specified.
if ~isempty(logr0)
    y0(2:3) = logr0;
end
if ~isempty(Sigma)
    y0(4:5) = sqrt([Sigma(1,1), Sigma(2,2)]);
    y0(6) = asin(Sigma(1,2)/sqrt(Sigma(1,1) * Sigma(2,2)));
end

opts.Display = 'off';
[y1,~,resid,~,~,~,J] = lsqnonlin(@(y) fun_pha(y)-x, y0, ...
    [0, -10, -10, 0, 0, -pi], [inf, 10, 10, 10, 3, pi], ... % search bounds
    opts); 

ci = nlparci(y1, resid, 'jacobian', J) - y1';
ci = 2 .* ci(:,2);  % 95% confidence interval

mu = [y1(2),y1(3)];
sigma = [y1(4),y1(5)];
Sigma = corr2cov(sigma,[1,sin(y1(6));sin(y1(6)),1]);

pha = Phantom('standard', grid, mu, Sigma);
pha.type = 'standard-fit';

N = y1(1); % scaling parameter denoting total number of particles
y_out = y1;


%-- Display mu, Sigma -----%
disp(' ');
disp('  p = ');
disp(pha.p);
%--------------------------%


tools.textheader();

end
