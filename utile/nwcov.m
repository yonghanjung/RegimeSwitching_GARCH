function V=nwcov(data,lag,meansub)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% covariance estimation using Newey-West
%  
% USAGE:
%   [V] = covnw(DATA)
%   [V] = covnw(DATA,NLAG,DEMEAN)
%
% INPUTS:
%   data   - vector of data 
%   lag   - the number of lag to use
%   mean - 0 or 1 if the mean should be subtracted when
%              computing the covariance 
%
% OUTPUTS:
%   V      - A K by K covariance matrix estimated 
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 5/1/2007

% Input Checking

T=size(data,1);
if nargin==1
    lag=min(floor(1.2*T^(1/3)),T);
    meansub=true;
elseif nargin==2
    meansub=true;    
end    
if isempty(lag)
    lag=min(floor(1.2*T^(1/3)),T);
end
if isempty(meansub)
    mean=true;
end
if ~ismember(meansub,[0 1]) 
    error('mean has to be a logical.')
end
if floor(nlag)~=nlag || nlag<0 
    error('lag has to be a integer.')
end
if ndims(data)>2
    error('error on data')
end

%%%%%%%%%%%%%%%%%%%%%%
if meansub
    data=data-repmat(mean(data),T,1);
end
 
w=(lag+1-(0:lag))./(lag+1);
V=data'*data/T;
for i=1:nlag
    Gammai=(data((i+1):T,:)'*data(1:T-i,:))/T;
    GplusGprime=Gammai+Gammai';
    V=V+w(i+1)*GplusGprime;
end


