function [L, loglik,pr,pt,h] = swgarchlikcoreH(parameters,P,data,k,ORDERS)
%
% Purpose:
%   use:  [L, loglik,pr,pt,h] = swgarchlikcoreH(parameters,P,data,k,ORDERS)
%   This function returns the likelihood for Haas Markov Switching GARCH
%   model. It also return the sum of indivudual likelihood and probability
%   state vector. 
%                                                                 
%   INPUT: - parameters: a vector of parameters. We want to evaluate the likelihood
%   function in these paramaters. The size of parameters depends on k, the number of regime
%              - P: transtion probability matrix. We want to evaluate the
%              likelihood with these probabilities of transition. Size of P
%              is k*k
%              - data: the data to compute the likelihood function
%              - k: the number of regimes
%              - ORDERS: orders of GARCH in the regime 
%
%   OUTPUT: - L: sum of the likelihood
%                 - loglik: vector of individual likelihood  
%                 - pr: vector of probabilities inference 
%                 - pt: vector of probabilities inference ex-ante
%
%  Author: Thomas CHUFFART
%  Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = ORDERS(2);
c = ORDERS(1);
q = ORDERS(3); 
m = max(p,q);
a = zeros(k,1);
b = zeros(k,1);

if c == 1,
    mu = parameters(1);
    epsi = data-mu;
else
    epsi = data;
end
T = length(epsi);
a0 = parameters(1+c:p+q+1:end);
for i = 1:p,
    a(:,i) = parameters(1+c+i:p+q+1:end);
end
for i = 1:q,
    b(:,i) = parameters(1+c+i+p:p+q+1:end);
end

b = diag(b);


A = [eye(k)-P;ones(1,k)];
I3 = eye(k+1);
c3 = I3(:,k+1);
if det(A'*A) ~= 0,
    pinf = ((A'*A)\A'*c3);
else
    pinf = ((A'*A + eye(k) * 1e-15)\A'*c3);
end

% fprintf('C3, pinf \n')
% P
% A
% c3
% pinf
% fprintf('\n')


pt = zeros(k,T+1);
pr = zeros(k,T);
h = zeros(k,T);
loglik = zeros(T,1);

h(:,1:m) = repmat(var(epsi),k,1);
pt(:,1) = pinf;

% fprintf('h,pt \n')
% h
% pt
% temp = normpdf(epsi(1,1),0,sqrt(h(1:k,1)));
% temp 
% fprintf('\n')

f(:,1)= normpdf(epsi(1,1),0,sqrt(h(1:k,1)));
LL(1:k,1)=pt(1:k,1).*f(:,1);
pr(:,1) = LL/sum(LL);
loglik(1,1)=log(ones(1,k)*LL);
pt(:,2) = P*pr(:,1);   

for t = m+1:T
   h(:,t) = a0+a*epsi(t-(1:p))^2+b*h(:,t-(1:q));
   LL(1:k,1)=pt(1:k,t).*normpdf(epsi(t,1),0,sqrt(h(1:k,t)));
   pr(:,t) = LL/sum(LL);
   pt(:,t+1)=P*pr(:,t);
   loglik(t,1)=log(ones(1,k)*LL);
end

pr = pr';
pt = pt';
L = -ones(1,T-1)*loglik(2:T);
loglik = ones(1,T-1)*loglik(2:T);