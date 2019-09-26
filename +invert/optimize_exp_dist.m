
% OPTIMIZE_EXP_DIST  Finds optimal lambda for exponential distance solver.
%=========================================================================%

function [x,lambda,out] = optimize_exp_dist(A,b,d_vec,m_vec,span,x_ex,Lex,x0,solver)
%-------------------------------------------------------------------------%
% Inputs:
%   A       Model matrix
%   b       Data
%   n       Length of first dimension of solution
%   lambda  Regularization parameter
%   span    Range for 1/Sf, two entry vector
%   x_ex    Exact distribution project to current basis
%   x0      Initial guess for solver    (Optional, default is zeros)
%   solver  Solver                      (Optional, default is interior-point)
%
% Outputs:
%   x       Regularized estimate
%   D       Inverse operator (x = D*[b;0])
%   Lx      Tikhonov matrix
%-------------------------------------------------------------------------%

%-- Parse inputs ---------------------------------------------------------%
if ~exist('x0','var'); x0 = []; end
if ~exist('x_ex','var'); x_ex = []; end
if ~exist('solver','var'); solver = []; end
%-------------------------------------------------------------------------%

lambda = logspace(log10(span(1)),log10(span(2)),70);

disp('Optimizing exponential distance regularization:');
tools.textbar(0);
for ii=length(lambda):-1:1
    out(ii).lambda = lambda(ii);
    [out(ii).x,~,out(ii).Lx,out(ii).Gpo_inv] = invert.exp_dist(...
        A,b,d_vec,m_vec,lambda,Lex,x0,solver);
    if ~isempty(x_ex); out(ii).chi = norm(out(ii).x-x_ex); end
    out(ii).Axb = norm(A*out(ii).x-b);
    
    tools.textbar((length(lambda)-ii+1)/length(lambda));
end

if ~isempty(x_ex)
    [~,ind_min] = min([out.chi]);
else
    ind_min = [];
end
lambda = out(ind_min).lambda;
x = out(ind_min).x;

end

