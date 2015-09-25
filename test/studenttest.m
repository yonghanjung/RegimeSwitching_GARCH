function results = studenttest(thetahat,data,stderr,test, alpha) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:
%   Perform a student test and return the t-sat and if we reject or not the
%   null hypothesis. 
%   We test H0: thetahat(k) = test
%   H1: thetahat(k) != test
%
% INPUT:
%   - thetahat: a vector of estimate coefficients
%   - data: The data require for the estimation
%   - stderr: The estimate standard errors
%   - test: The value wich is used for the test, if no input, test = 0 
%   - alpha: H0 is accepted if |tk| < t0,005 with probability alpha, by
%     default, alpha = 5%
%
% OUTPUT:
%   -results, a structure wich contains:
%       tstat: The t-satitistic
%       pval: P-value 
%       h: 0 if H0 is accepted, 1 if H0 is rejected
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin' INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,
    error('wrong input argument');
end

if nargin < 4,
    alpha = 0.05;
end

if isvector(thetahat) ~= 1
    error('theta has to be a vector');
end



if isvector(stderr) ~= 1
    error('theta has to be a vector');
end

if length(thetahat) ~= length(stderr),
    error('thetahat et stderr have to be equals');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = length(thetahat);
tstat(1:t) = (thetahat(1:t)-test)./stderr(1:t);

dim = length(data);
pval = zeros(t,1);
h = zeros(t,1);

for i = 1:t,
    pval(i) = 1-tcdf(abs(tstat(i)),dim-t);

    if pval(i) < alpha,
        h(i) = 0;
    else
        h(i) = 1;
    end
end

results.tstat = tstat;
results.pval = pval; 
results.H = h; 
end
    