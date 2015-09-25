function [sumLLF LLF h epsihat] = argarchLik(para,y,ORDERS,others,nbcoef,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function computes the likelihood of an AR(P)-GARCH(p,q) process.
%   Depend on the type of GARCH model you choose. Can be GARCH(p,q),
%   LST-GARCH(p,q) of GJR-GARCH(p,q,o)
%
% INPUTS:
%   para: vector of the parameters
%   y: vector of the data 
%   ORDERS: A vector of length = 4 with ORDERS(1) = C (=1 if constante in
%                 the AR process), ORDERS(2) = P (AR order), ORDER(3) = p (ARCH order)
%                 and ORDER(4) = q (GARCH order)
%   others: if type = 2, others is the smoothing function
%              if type = 3 others is the order of asymetrie 
%   nbcoef: 
%   type: Type of garch component: 1: GARCH, 2:LST-GARCH, 3: GJR-GARCH. More soon
%
% OUTPUTS:
%   sumLLF: The sum of the individual likelihood
%   LLF: The likelihood in theta 
%   h: h estimated
%   epsihat: esimated residuals
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 6,
    error('wrong input argument');
end

if size(para,1) < size(para,2)
    para = para';
end

if length(ORDERS) ~= 5,
    error('ORDERS has to be a vector of length 5');
end

if type ~= 4, 
    c = ORDERS(1);
    P = ORDERS(2);
    Q = ORDERS(3);
    p = ORDERS(4);
    q = ORDERS(5);
    if type == 3,
        assym = others;
    else
        assym = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = length(y);
nbcoefar = Q+P+c; 
maxi = max(ORDERS);

if c == 0 && P == 0 && Q == 0,
    nbparamu = 0;
    paramu = 0;
else
    paramu = para(1:P+c+Q);
    nbparamu = length(paramu);
end

paraGARCH = para(nbparamu + 1: end); 

if type == 1 || type == 3,
    paraGARCH = para(nbcoefar + 1: end); 
elseif type == 2,
     paraLSTGARCH = para(nbcoefar + 1: end-1);
     theta = para(nbcoef);
end

h = zeros(length(y),1);
h(1:maxi) = var(y);
mu = zeros(dim,1);

if type == 1,
    for t = (maxi + 1):dim,
        if c == 0 && Q == 0 && P ==0,
            mu(t) = 0;
        else
            mu(t) = paramu'*[1:c; y(t-(1:P)); y(t-(1:Q))-mu(t-(1:Q))];
        end
        h(t) = paraGARCH' * [1 ; (y(t-(1:p)) - mu(t-(1:p))).^2 ; h(t-(1:q))];
    end
elseif type == 2
    epsi = zeros(dim,1);
    F = zeros(dim,1);
    if strcmp( 'L', others),
        for t = (maxi + 1):dim,  
            if paramu == 0,
                mu(t) = 0;
            else
                mu(t) = paramu'*[1:c; y(t-(1:P)); y(t-(1:Q))-mu(t-(1:Q))];
            end
            epsi(t) = y(t) - mu(t);
            F(t) = (1/(1+exp(-theta*(epsi(t-(1:p)))))-0.5);
            h(t) = paraLSTGARCH' * [1; epsi(t-(1:p)).^2 ;F(t)*(epsi(t-(1:p))).^2 ; h(t-(1:q))];
        end
    else
        invtheta = 1/theta; 
        for t = (maxi + 1):dim,
            mu(t) = paramu'*[1:c; y(t-(1:P)); y(t-(1:Q))-mu(t-(1:Q))];
            h(t) = paraLSTGARCH' * [1 ; (y(t-(1:p)) - mu(t-(1:p))).^2 ; expcdf((y(t-(1:p)) - mu(t-(1:p))),invtheta)*(y(t-(1:p)) - mu(t-(1:p))).^2  ; h(t-(1:q))];
        end
    end
elseif type == 3,
    for t = (maxi+1):dim,
        if paramu == 0,
            mu(t) = 0;
        else
            mu(t) = paramu'*[c; y(t-(1:P)); y(t-(1:Q))-mu(t-(1:Q))];
        end
            h(t) = paraGARCH' * [1 ; (y(t-(1:p)) - mu(t-(1:p))).^2  ; h(t-(1:q));(y(t-(1:assym))-mu(t-(1:assym))).*(y(t-(1:assym))-mu(t-(1:assym))).^2];
    end
end

t = 1:dim;
LLF = 0.5*((log(h(t))) + (((y(t)-mu(t)).^2)./h(t)) + log(2*pi));
sumLLF = (sum(LLF)); 
if isnan(LLF)
    sumLLF=1e+6;
end

h = h(t,1) ;
epsihat = (y(t)-mu(t))'; 

end