function [thetahat results struct]= swgarchest(data,flag,ORDERS,reg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function estimates a Markov-Switching GARCH(p,q) process.
%   You can choose the type of the MS-GARCH process: Haas & all (2004) or
%   Klaassen (2002).
%
% INPUTS:
%   data: The data 
%   flag: type of model to estimate. 1: Klaassen (2002) 2: Haas (2004) 
%   %%% ADDED %%% ORDERS should be [1,1,1] for now %%% ADDED %%%
%   ORDERS: A vector of length = 3 with ORDERS(1) = c the conditional mean, 
%   ORDERS(2) = p (GARCH order), ORDER(3) = q (GARCH order)

    % ORDER = [1,1,1] or [0,1,1] 

%   reg: the number of regime

    % reg = 2 
    % nargin = 4 

%   startval: vector of initial values for the likelihood initialization 
%   startM: transition probability matrix starting values for the likelhood
%   initialization
%   hreal: The real value of h if they are known or a unbiaised estimator
%   like the realized volatility
%
% OUTPUTS:
%   thetahat: Vector of estimate coefficient
%   results: a structure with many statistics (AIC, BIC, loss functions...)
%   struct: a structure wich contains output of fmincon
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check the inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 7,% This is My case! 
    flagmse = 0;
else
    flagmse = 1;
    if size(hreal,2) > 1,
        error('hreal vector should be a column vector');
    end
end
    
if nargin < 4 
    error('number of regime mispecified');
elseif length(reg) > 1 || reg < 0
    error('reg has to be a positive scalar');
end

k = reg; 

if nargin < 3,
    error('ORDERS is missing');
elseif length(ORDERS) ~= 3,
    error('ORDERS has to be a vector of length 3');
else
    c = ORDERS(1);
    p = ORDERS(2);
    q = ORDERS(3);
    if (length(p) > 1) || p < 0 || length(q) > 1|| q < 0
        error('p and q should be positive scalars')
    end
end

if nargin < 2, 
    flag = 2;
end

if nargin < 1,
    error('you should read the abstract')
end

% This is My case! 
if nargin < 6 
%     startM = diag(ones(k));
    startM = diag(ones(k,1)); % Identity matrix as the first 
elseif size(startM,1) ~= k || size(startM,2) ~= k
    error('startM has to be a k*k matrix')
end

% Check whether the column sum of startM is equal to 1 
A = sum(startM);
if A ~= ones(1,k),
    error('startM has to be a probability transtion matrix')
end



% MY CASE! 
if nargin < 5 || isempty(startval),
    startval = 0; % What the hell is that?
end


% MY CASE! 
if size(startval,1) ~= k || size(startval,2) ~= (p+q+c+1), 
    if startval == 0,
        % Printed
    fprintf('\n Valdep is empty, we assign it default value \n'); 
    
    else
        fprintf('\n startval is mispecified, we assign it default values \n')
    end
    
    % if [1,1,1] 
    if c ~= 0,
            [Y X] = matrixlag(data,0,c);
            bhat = X\Y; 
            bhat = real(bhat);
    
    % if [0,1,1] 
    else
        bhat = [] ;
    end
    
    alpha  =  .05*ones(k,p);
    beta   =  0.05*ones(k,q)/q;
    omega  =  0.01*ones(k,1);  
    
    % if [0,1,1] 
    if isempty(bhat),
        startval = [omega ; alpha ; beta];
    
    % if [1,1,1] 
    else
        if p ~= 0,
            startval = [bhat ; omega ; alpha ; beta]; % num: 7
        else
            startval = bhat;
        end
    end
end

if size(startval,2)>1
    startval = startval';
    
end

if size(data,2) > 1 
   error('Data vector should be a vector');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = ORDERS(1);
p = ORDERS(2);
q = ORDERS(3);
k = reg;
% Number of GARCH parameters, 2*(1+1+1 (omega)  ) = 6 
nbparaGarch = k*(p+q+1); % =6 

% Number of paramters, 8 
nbpara = k*(p + q + 1 + c);

if k == 2,

%     startvaltot = [reshape(startval,nbpara,1) ; startM(1,1);startM(2,2)];

% 9 dim if [1,1,1]
% [bhat; omega; alpha; beta; startM(1,1); startM(2,2)] 
    startvaltot = [reshape(startval,length(startval),1) ; startM(1,1);startM(2,2)];
else
    startvaltot = [reshape(startval,nbpara,1) ; reshape(startM,k^2,1)];
end

% Lower Bound: [0,0,0,0,0,0,0,0,0] 
LB = zeros(1,length(startvaltot));

% No Upper bound 
UB = [];
y = data;

options  =  optimset('fmincon');
options  =  optimset(options , 'Algorithm ','interior-point');
% options  =  optimset(options , 'Algorithm ','active-set');
options  =  optimset(options, 'Hessian','bfgs');
options  =  optimset(options , 'TolFun'      , 1e-006);
options  =  optimset(options , 'TolX'        , 1e-006);
options  =  optimset(options , 'TolCon'      , 1e-006);
options  =  optimset(options , 'Display'     , 'iter');
options  =  optimset(options , 'Diagnostics' , 'on');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxIter'     , 1500);
options  =  optimset(options , 'Jacobian'     ,'off');
options  =  optimset(options , 'MeritFunction'     ,'multiobj');
options  =  optimset(options , 'MaxFunEvals' , 3000);

if k == 2,
%     a = -eye(8);
    % a = -I9
    % for all param, 0 <= param <= 1 for stationarity 
    a = -eye(length(startvaltot));
    
    % A = [-I9; +I9], 18 by 9 
    A = [a;-a];
%     bt = [zeros(1,8) ones(1,8)]
    bt = [zeros(1,length(startvaltot)) ones(1,length(startvaltot))]; 
    
    % bt = [0;....;0; 1;....;1] 
    Aeq = [];
    beq = [];
else
    a = -eye(length(startvaltot));
    A = [a;-a];
    bt = [zeros(1,length(startvaltot)) ones(1,length(startvaltot))] ;
    Aeq = [zeros(k,nbparaGarch) blokdiag(ones(1,k),k)];
    beq = ones(1,k);
end

% fprintf('A,bt')

% Convex optimization 
% Maximize the likelihood function 
% Control variable is x 

% swgarchlik(x,data,reg=2,[0,1,1],2),x0(startvalot),A,bt,Aeq,beq,LB,UB,
% nonlinearconstraint : constMSGARCH 
[thetahat,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) swgarchlik(x,data,reg,ORDERS,flag),startvaltot,A,bt,Aeq,beq,LB,UB,@(x) constrMSGARCH(x,k,nbpara),options); 

% [thetahat,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) swgarchlik(x,data,reg,ORDERS,flag),startvaltot,[],[],[],[],[],[],@(x) constrMSGARCH(x,k,nbpara),options); 
[LLF,likelihoods,~,p,pt,smoothprob,h] = swgarchlik(thetahat,data,reg,ORDERS,flag);

H = zeros(length(h),1);
for i = 1:length(h),    
    H(i,1) = p(i,:)*h(:,i);
end


struct.output = output;
struct.fval = fval;
struct.exit = exitflag;
struct.lamb = lambda;
struct.grad = grad;
struct.hessian= hessian;


[VCV]=robustvcv(@swgarchlik,thetahat,0,y,reg,ORDERS,flag);
stderr=sqrt(diag(VCV));

results.stderr = stderr;
results.LLF = LLF; 
results.lik = likelihoods;
results.p = p;
results.pt = pt;
results.smooth = smoothprob;
results.AIC = InfoCrit(-LLF,thetahat,data,1);
results.BIC = InfoCrit(-LLF,thetahat,data,2);
results.H = H;
results.h = h;

if flagmse,
    dim = length(data);
    results.mise = lossfun(hreal(1:dim),H,2);
    results.qlike = lossfun(hreal(1:dim),H,3);
    results.mae = lossfun(hreal(1:dim),H,4);
    
    results.mise = lossfun(data(1:dim).^2,H,2);
    results.qlike = lossfun(data(1:dim).^2,H,3);
    results.mae = lossfun(data(1:dim).^2,H,4);
end

end

