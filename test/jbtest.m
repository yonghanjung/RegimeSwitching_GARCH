function results = jbtest(data, alpha, k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function compute the Jarque and Bera's test. 
% This test is used for test the normality of the data. 
% H0: The data are normally distributed
% H1: The data are not normally distributed

% In fact, the Jarque and Bera's test focus on the kurtosis and the
% skewness of the data. Data which are normally distributed have an
% expected skewness equals to 0 and an expected kurtosis equals to 3. 
% So the null hypothesis is:
% H0: S=0 and K=3
%We will accept this hypothesis if the p-value is less than alpha.
%
%   INPUT:
%       - data: The data that we want test
%       - alpha: The level of the test
%       - k: Number of explanatory variables if the data come from the
%             residues of a linear regression. Otherwise, k = 0.
%
%   OUTPUT:
%       results, a structure wich contain:
%           - stat: The JB's stat
%           - pval: The p-value
%           - H: 0 if H0 is accepted, 1 if H0 is rejected
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin' INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,
    alpha = 0.05;
end

if nargin<3,
    k = 2;
end

if isvector(data) ~= 1,
    error('data has to be a vector');
end

if size(data,1)<size(data,2)
    data = data';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = length(data);

S2 = skewness(data)^2/6;
K = kurtosis(data); 

B = ((K-3)^2)/24;

stat = (dim-k)*(S2+B); 
pval = 1-chi2cdf(stat,2);
H = pval<alpha;

results.stat = stat;
results.pval = pval;
results.H = H;

end

