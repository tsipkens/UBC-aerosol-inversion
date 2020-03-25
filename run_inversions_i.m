
% RUN_INVERSIONS_I  Optimize exponential distance regularization w.r.t. lambda.
% Author: Timothy Sipkens, 2019-05-28
%=========================================================================%

if iscell(phantom.Sigma)
    Gd = phantom.Sigma{1};
else
    Gd = phantom.Sigma;
end

%-- Gd properties -----------------%
l1 = sqrt(Gd(1,1));
l2 = sqrt(Gd(2,2));
R12 = Gd(1,2)/(l1*l2);
Dm = Gd(1,2)/Gd(2,2); % s1*R12/s2
%----------------------------------%

[x_ed_lam,lambda_ed_lam,out_ed_lam] = ...
    optimize.exp_dist_op(...
    Lb*A,Lb*b,[0.1,10],Gd,...
    grid_x.elements(:,2),grid_x.elements(:,1),x0);
disp('Process complete.');
disp(' ');

err.ed = norm(x_ed_lam-x0);


