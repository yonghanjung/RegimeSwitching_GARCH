function [processroots, absroots] = arroots(phi, p, c)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function calculates the roots of an AR(P) process.
%    
% INPUTS:
%   phi: AR parameters, a vector of (p+c)*1
%   p: AR process order
%   c: constant order
%   
% OUTPUTS:
%   processroots: the roots of the AR process
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUTS CHECKIN' %%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 3,
    error('Wrong input argument');
end

if size(phi,2)> size(phi,1)
    phi = phi';
end
if size(phi,2)~=1
    error('Phi is not a column vector');
end

if size(phi,1)~= p+c,
    error('Lenght of phi not consistant with the order of the AR')
end

if isvector(p)~=1 || isscalar(c)~=1,
    error('p or c have to be scalar');
end

if p <0 ||  c <0,
    error('p or c have to be postive');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute AR-ROOTS

if c ~=0,
    arphi = phi(c+1:p+c);
else
    arphi = phi(1:p); 
end

ar([1 p]) = arphi;
processroots=roots([-ar 1]);
absroots=abs(processroots);

end





    
    