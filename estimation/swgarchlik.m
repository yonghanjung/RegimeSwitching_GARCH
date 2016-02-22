function [LLF,likelihoods,errors,p,pt,smoothprob,h] = swgarchlik(para,data,k,ORDERS,flag)

% Function to calculate the loglikelihood for different
% Markov Switching Regime GARCH (MRS-GARCH) model
%
% INPUT:
% para:           vector of parameters, ordered by regime. The last parameters
%                   are the transition probabilities
% data:           vector of data
% k                 The number of regime
% ORDERS:     The order of the regime switching garch models. For the
%                   moment, all the GARCH models has the same
%                   specification, we will include soon the possibility to have
%                   a GARCH(1,1) and a GARCH(1,2) for exemple
%flag:             = 1 Klaassen's (2002, Empirical Economics) specification
%                   = 2 Haas & all (2004, Journal of Financial Econometric)
%                   specification
%                  
%
% OUTPUT:
% LLF:               value of log-likelihood function
% likelihoods:     vector of log-likelihoods
% errors:            vector of errors
% p:                  Matrix of filtered probabilities, the colum k is the vector of probabilities 
%                       to be in regime k at time t given the information a
%                       time the information at time t
% pt:                 Matrix of ex-ante probabilities, the colum k is the vector of probabilities 
%                       of being in regime k at time t given the information
%                       at time t-1
% smoothprob:   Matrix of smoothed probabilities. The column k is the vector of probabilitie 
%                       of being in regime k at time t given all the
%                       information in the sample
% h:                   Estimated volatility 
%
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization of the loglikelihood function and the other variables

p = ORDERS(2);
c = ORDERS(1);
q = ORDERS(3); 
% nbparatot = (p+q+c+1)*k;
nbparatot = (p+q+c)*k ; % (1+1+1) * 2 = 6 

% [hb]

% [bhat; omega; alpha;; beta] 
% prob = param(7,8)
% parameters = [omega, alpha, beta] 
% parameters = (1,2,3,4,5,6)
prob = para(nbparatot+1:end); 
parameters = para(1:nbparatot);



% fprintf('prob \n')
% prob
% parameters
% fprintf('\n')

if k > 2,
    P = reshape(prob,k,k);
else
    P = [prob(1) 1-prob(2) ; 1-prob(1) prob(2)];
end
T = length(data);

% We use the m function swgarchlikcore

switch flag
    case 1 % Klaassen's specification
        [L,loglik,p,pt,h] = swgarchlikcoreK(parameters,P,data,k,ORDERS);         
        L
    case 2 %Haas specification
        [L,loglik,p,pt,h] = swgarchlikcoreH(parameters,P,data,k,ORDERS);            
end

% Smoother to calculate the smoothing probabilities
% KIM's Smoother for smooth probabilities 

Xitp1t = pt;
Xitt = p;
smoothprob = Xitt;
for i = 1:(T-1)
    smoothprob(T-i,:) = ((P'*(smoothprob(T-i+1,:)'./(Xitp1t(T-i+1,:)'))).*Xitt(T-i,:)')';
end

LLF = L;
likelihoods = loglik;
errors = 0;
end