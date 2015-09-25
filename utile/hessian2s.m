function H = hessian2s(f,arg,varargin)
%
% INPUTS:
%   f         - Function name, fval = func(x,varargin)
%   arg            - Vector of parameters (N x 1)
%   VARARGIN     - Optional arguments passed to the function
%
% OUTPUTS:
%   H            - Finite differnce, 2-sided hessian
% COMMENTS:
%   Modification of hessian() from the jpl toolbox


% Code originally from COMPECON toolbox [www4.ncsu.edu/~pfackler]
% documentation modified to fit the format of the Econometrics Toolbox
% by James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jlesage@spatial-econometrics.com
%
% Further modified (to do 2-sided numerical derivs, rather than 1) by:
% Modifications Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 3    Date: 9/1/2005
%
% I have make some modificiations which seems to not change results
% Modifications Copyright 2: Thomas Chuffart 


T = size(arg,1);

if size(arg,2)~=1
    error('X must be a column vector.')
end

try
    feval(f,arg,varargin{:});
catch FE
    error(['There was an error evaluating the function:' FE.message]);
end

fx = feval(f,arg,varargin{:});

h = eps.^(1/3)*max(abs(arg),1e-2);
argh = arg+h;
h = argh-arg;
ee = sparse(1:T,1:T,h,T,T);

gp = zeros(T,1);
gm = zeros(T,1);

for i=1:T
    gp(i) = feval(f,arg+ee(:,i),varargin{:});
    gm(i) = feval(f,arg-ee(:,i),varargin{:});
end

hh=h*h';
Hm=NaN*ones(T);
Hp=NaN*ones(T);
H = zeros(T);

for i=1:T
    for j=i:T
        Hp(i,j) = feval(f,arg+ee(:,i)+ee(:,j),varargin{:});
        Hp(j,i)=Hp(i,j);
        Hm(i,j) = feval(f,arg-ee(:,i)-ee(:,j),varargin{:});
        Hm(j,i)=Hm(i,j);
        H(i,j) = (Hp(i,j)-gp(i)-gp(j)+fx+fx-gm(i)-gm(j)+Hm(i,j))/hh(i,j)/2;
        H(j,i) = H(i,j);
    end
end

