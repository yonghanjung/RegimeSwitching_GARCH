function [y,x]=matrixlag(x,nb_lags,const)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the lag of a vector nb_lags times
%
%
%   INPUT:
%       - data: x, the data 
%       - nb_lags: The number of lags 
%       - const: 1 if we want a constant in the lag of x 
%
%   OUTPUT:
%       [X Y] a matrix with
%           - X the vector x 
%           - Y: A matrix n * nb_lags 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin' INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,
    error('Wrong inputs arguments');
end

if nargin < 3,
    const = 1;
end

if isscalar(nb_lags ) ~= 1,
    error('nb_lags has to be a scalar');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_lags=nb_lags+1;
dim=size(x,1);

X=[x;zeros(nb_lags,1)];
lagmatrix=repmat(X,nb_lags,1);
lagmatrix=reshape(lagmatrix(1:size(lagmatrix,1)-nb_lags),dim+nb_lags-1,nb_lags);
lagmatrix=lagmatrix(nb_lags:dim,:);
% size(lagmatrix)
y=lagmatrix(:,1);
x=lagmatrix(:,2:nb_lags)
if const==1
    x=[ones(size(x,1),1) x];
end


