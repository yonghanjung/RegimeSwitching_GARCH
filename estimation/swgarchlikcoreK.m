function [L, loglik,pr,pt,h] = swgarchlikcoreK(parameters,P,data,k,ORDERS)
%
% Purpose:
%   use:  [L, loglik,pr,pt,h] = swgarchlikcoreK(parameters,P,data,k,ORDERS)
%   This function returns the likelihood for Klaassen Markov Switching GARCH
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
%e
%   OUTPUT: - L: sum of the likelihood
%                 - loglik: vector of individual likelihood  
%                 - pr: vector of probabilities inference 
%                 - pt: vector of probabilities inference ex-ante
%                 - h: the true volatility
%
%  Author: Thomas CHUFFART
%  Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = ORDERS(2);
c = ORDERS(1);
q = ORDERS(3); 
m = max(p,q);
T = length(data);
a0 = parameters(1+c:p+q+1:end);
a = zeros(k,1);
b = zeros(k,1);

if c == 1,
    mu = parameters(1);
    epsi = data-mu;
else
    epsi = data;
end

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
pinf = ((A'*A)\A'*c3);

pt = zeros(k,T+1);
pr = zeros(k,T);
h = zeros(k,T);
loglik = zeros(T,1);
pij = zeros(k);

h(:,1:m) = repmat(var(epsi),k,1);
pt(:,1) = pinf;
f(:,1)= normpdf(epsi(1,1),0,sqrt(h(1:k,1)));
LL(1:k,1)=pt(1:k,1).*f(:,1);
pr(:,1) = LL/sum(LL);
loglik(1,1)=log(ones(1,k)*LL);
pt(:,2) = P*pr(:,1);   
for i = 1:k,
    for j = 1:k
        pij(i,j) = P(i,j)*pr(j,1)/pt(i,2);
    end
end  
hit1 = pij * h(:,1);
v(:,1) = hit1;

for t = m+1:T
   h(:,t) = a0+a*epsi(t-(1:p))^2+b*v(:,t-(1:q));
   f(:,t)= normpdf(epsi(t,1),0,sqrt(h(1:k,t)));
   LL(1:k,1)=pt(1:k,t).*f(:,t);
   pr(:,t) = LL/sum(LL);
   pt(:,t+1)=P*pr(:,t);
   loglik(t,1)=log(ones(1,k)*LL);
    for i = 1:k,
        for j = 1:k
            pij(i,j) = P(i,j)*pr(j,t)/pt(i,t+1);
        end
    end  
   hit1 = pij * h(:,t);
   v(:,t) = hit1;
end


pr = pr';
pt = pt';
L = -ones(1,T)*loglik(1:T);
loglik = ones(1,T)*loglik(1:T);