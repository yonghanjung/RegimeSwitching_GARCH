function results = ljungboxtest(res,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:
%   This function tests the autocorellation of the squared residuals. It
%   computes the Q-statistic discover by Ljung and Box in 1978.
%
%   G. M. Ljung; G. E. P. Box (1978). "On a Measure of a Lack of Fit in Time
%   Series Models". Biometrika 65 (2): 297?303
%
% 
% INPUTS:
%   - res: A vector of residuals
%   - p: The number of lag we want to test autocorrelation
%   - alpha: The level of the test, is alpha is not defined alpha = 5%
%
% OUTPUTS:
%   - results, a structure wich contains:
%       - stat: The Q-stat
%       - pval: The p-value of the test
%
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin' INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,
    error('2 inputs are needed');
end

if isvector(res)~=1,
    error('res has to be a vector');
end

if size(res,1)<size(res,2)
    res = res';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stat = zeros(p,1);
dim = length(res);

rhoemp = autocorr(res,20,0);

for k=1:p
    stat(k)=dim*(dim+2)*sum(rhoemp(1:k).^2./(dim-(1:k)'));
end

pval = 1-chi2cdf(stat,(1:k)'); 
results.pval = pval;
results.stat = stat; 

end
