function results = testEngel(x,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:
%   This function computes the Engel test. It test the presence of ARCH
%   process in data's residuals.
%
%   Consider the regression:
%   epsi^2(t) = w + a1*epsi^2(t-1) ... + a2*epsi^2(t-p)
%   (epsi can to come from an estimation of an AR process...)
%
%   It follow two steps:
%       - Step 1: Estimate the regression below
%       - Step 2: Compute the R2 and the test statistic.
%
% INPUTS:
%       - x : The data wich are supected to follow an ARCH process
%       - p: The number of lag
%
% OUPUTS:
%       results, a structure wich contain:
%           - stat: The statistic of the test = T*R2
%           - pval: The p-value of the test
%           - H: 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checkin' INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin ~= 2)
   error('Wrong number of arguments to arch');
end;

x=x.^2;

[y,X] = matrixlag(x,p,1);
bhat = X\y; 

yhat = X * bhat;
epsilon = y - yhat;
ytilde = y - mean(y);
R2= 1 - (((epsilon')*(epsilon))/((ytilde')*(ytilde)));
stat = length(y)* R2;
pval = 1-chi2cdf(stat,p);
H = pval < 0.05;

results.stat = stat;
results.pval = pval;
results.H = H;

end
