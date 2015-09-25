function gama = autocov(X,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function computues empirical auto-covariance of the data
%
% INPUTS:
%   - X: The data
%   - k: The number of lags
%
% OUTPUTS:
%   - gama: The empirical auto-covariance 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin~=2
    error('2 inputs needed')
end

T = length(X);
gama = mean(X(1:T-k).*X(k+1:T)) - mean(X)*mean(X(k+1:T));

end 