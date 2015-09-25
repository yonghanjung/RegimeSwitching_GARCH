function [VCV,A,B,scores,hess]=robustvcv(fun,theta,nw,varargin)
% Compute Robust Variance Covariance matrix numerically
%
% USAGE:
%     [VCV,A,B,SCORES,HESS,GROSS_SCORES]=robustvcv(FUN,THETA,NW,VARARGIN)
%
% INPUTS:
%     FUN           - Function name ('fun') or function handle (@fun) which will
%                       return the sum of the log-likelihood (scalar) as the 1st output and the individual
%                       log likelihoods (T by 1 vector) as the second output.
%     THETA         - Parameter estimates 
%     NW            - Number of lags 
%     VARARGIN      - Other inputs to the likelihood function
%
% OUTPUTS:
%     VCV           - Estimated robust covariance matrix (see White 1994)
%     A             - A portion of robust covariance
%     B             - B portion of robust covariance
%     SCORES        - T x num_parameters matrix of scores
%     HESS          - Estimated Hessian (Expectation of second derivative)
%
% COMMENTS:
%     This function simplifies calculating sandwich covariance estimators for (Q)MLE estimation

% Input Checking

if size(theta,1)<size(theta,2)
    theta=theta';
end

%Computation of the VCV matrix

k=length(theta);
h=diag(max(abs(theta*eps^(1/3)),1e-8));

[~,likelihood]=feval(fun,theta,varargin{:});

T=length(likelihood);

LLFp=zeros(k,1);
LLFm=zeros(k,1);
likelihoodp=zeros(T,k);
likelihoodm=zeros(T,k);

% Computation of the likelihood near the optimal solution
for i=1:k
    thetap=theta+h(:,i);
    [LLFp(i),likelihoodp(:,i)]=feval(fun,thetap,varargin{:});
    thetam=theta-h(:,i);
    [LLFm(i),likelihoodm(:,i)]=feval(fun,thetam,varargin{:});
end
% Computation of first derivative 
scores=zeros(T,k);
h = diag(h);
for i = 1:k,
    scores(:,i)=(likelihoodp(:,i)-likelihoodm(:,i))./(2*h(i));
end
hess=hessian2s(fun,theta,varargin{:});

A=hess/T;
hess=A;
Ainv=A^(-1);
if nw==0
    B=cov(scores);
    VCV=(Ainv*B*Ainv)/T;
else
    B=nwcov(scores,nw);
    VCV=(Ainv*B*Ainv)/T;
end