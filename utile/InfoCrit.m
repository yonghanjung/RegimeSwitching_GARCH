function [crit1,crit2] = InfoCrit(LLF,para,data,flag)

% PURPOSE: compute information criterions
%
% INPUT:
%   LLF:            value of the log likelihood at the optimum
%   para:           vector of parameter estimate
%   flag             1 for AIC, 2 for BIC....
%                  
%
% OUTPUT:
%   crit:               value of the criterions
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3,
    error('wrong number of input');
end

[m,n] = size(para);
if n ~= 1 && m ~= 1, 
    error('para has to be a voector');
end

if (length(flag) > 1 || length(LLF) > 1) 
        error('flag and LLF should be scalars');
end

k = length(para);
n = length(data);
switch flag
    case 1
        crit1 = (2*k - 2*LLF);
        crit2 = (2*k - 2*LLF) + (2*k*(k+1))/(n-k-1);
    case 2    
        crit1 = k*log(n) - 2 * LLF;
end
