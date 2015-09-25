function [simulateEpsi,H,Y] = argarchSim(type,dist,dim,para,ORDERS,others)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function simulates an univariate AR-GARCH process with
%   Gaussian or student innovations. This function generates 1000 observations in more
%   for reducing any bias. You can choose the type of GARCH pocesses
%   between LST-GARCH, GARCH or GJR-GARCH for the moment. 
%
% INPUTS:
%   type: Type of garch component: 1: GARCH, 2:LST-GARCH, 3: GJR-GARCH. More soon
%   dist: Distribution of the error term 'GAUSSIAN' or 'T' (The first only
%   is available for the moment)
%   dim: length of the process simulated
%   para: The vector of the paramaters of the process ->
%                if type = 1 -> [const ARpara MApara ARCHpara GARCHpara]
%                if type = 2 -> [const ARpara MApara ARCH1para ARCH2para GARCHpara theta]
%                if type = 3 -> [const ARpara MApara ARCHpara GARCHpara GJRpara]
%   ORDERS: A vector of length = 5 with ORDERS(1) = C (=1 if constante in
%                 the AR process), ORDERS(2) = P (AR order), ORDERS(3) = Q (MA order) 
%                 ORDER(4) = p (ARCH order) and ORDER(5) = q (GARCH order)
%   others: if type = 2 equals to 'E' or 'L" -> Define the type of the smooth function
%               if type = 3, order of the asymetrie    
%
% OUTPUTS:
%   simulateEpsi: GARCH errors
%   H: Conditional variance
%   Y: Simulated process
%
% Author: Thomas Chuffart
% Mail: thomaschuffart@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Checkin' INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 5,
    error('Wrong input number');
end

if isscalar(dim),
    T = 1000+dim; 
else
    error('dim has to be a scalar');
end

if length(ORDERS) ~= 5,
    error('ORDERS has to be a vector of length 5');
else 
    c = ORDERS(1);
    P = ORDERS(2);
    Q = ORDERS(3);
    p = ORDERS(4);
    q = ORDERS(5);
end

if isscalar(type),
    if type > 4 || type <1,
    error('type has to be equal at 1, 2 , 3 or 4');
    end
else
    error('type has to be a scalar');
end

if isempty(q),
    q = 0;
end
    
if nargin <  6,
    others = 0 ;
end

if (length(p) > 1) || (length(q) > 1) || p < 0 ||q < 0
    error('p and q have to be scalar > 0');
end

if (length(P) > 1) || (length(Q)>1) || P < 0 || Q < 0
    error('P and Q have to be scalar');
end

if (length(c) > 1) || c < 0
    error('c has to be scalar');
end

if p == 0 ,
    q = 0; 
end 

if size(para,2)>size(para,1)
    para = para';
end

if type == 1,
    if length(para) ~= sum(ORDERS) +1,
                error(['Error in the number of parameters. You need a vector of length ' num2str(sum(ORDERS)+1)]);
    end
end

if type == 2,
    if ischar(others),
        EST = strcmp('E',others) ;
        LST = strcmp('L',others) ;
        if EST == 0  &&  LST == 0, 
            error('Erreur in smooth function input, has to be equal to a string variable L or E');
        end
    end
    if length(para) ~= (2*p)+P+q+Q+c+2,
        error(['Error in the number of parameters. You need a vector of length ' num2str((2*p)+P+q+Q+c+2)']);
    end
end

 if type == 3 
     if isscalar(others),
         if length(para) ~= sum(ORDERS)+1+others,
            error(['Error in the number of parameters. You need a vector of length ' num2str(sum(ORDERS)+1+others)']);
         else
             assym = others; 
             maxi = max(ORDERS);
         end
     else
         error('Error, others is not a scalar');
     end
 end
 
 if type == 1 || type == 2,
     maxi = max(ORDERS);
 end
 
 if strcmp(dist,'GAUSSIAN') 
    n = randn(T+maxi,1);
elseif strcmp(dist,'T')
    n = trnd(para(end), T);
else
    error('bad distribution');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1; 
if c == 0,
    a = 0;
else
    a = para(i);
    i = 1+i;
end

if P == 0,
    b = 0;
    P = 1;
else
    b = para(i:P+i-1);
    i = P + i; 
end

if Q == 0,
    d = 0;
    Q = 1;
else
    d = para(i:Q+i-1);
    i = Q + i; 
end

if p ~= 0,
    w = para(i);
    i = i +1; 

    if type == 1 || type == 3
        alpha = para(i:i+p-1);
        i = i + p ;
    elseif type == 2,
        alpha1 = para(i:i+p-1);
        i = i + p ;
        alpha2 = para(i:i+p-1);
        i = i + p ;
        theta = para(end);
    end

    if q == 0, 
        beta = 0; 
        q = 1; 
    else
        beta = para(i:i+q-1);
        i = i + q;
    end

    if type == 3, 
        gamma = para(i:end);
    end 

    if type == 1,
        para = [w ;alpha; beta];
        var = w / (1 - sum(alpha) - sum(beta)) ; 
    elseif type == 2 
        para = [w ; alpha1 ; alpha2 ; beta];
        var = w / (1 - sum(alpha1) - sum(beta)) ; 
        if var<0,
            var = 0.5;
        end
    elseif type == 3, 
        para = [w ; alpha; beta; gamma];
        var = w / (1 - sum(alpha) - sum(beta)-0.5*sum(gamma)) ; 
    end
else 
    para = 0;
end

h = 0.5 * ones(T+maxi,1);
epsi = sqrt(var) * ones(T+maxi,1); 
y = zeros(T+maxi,1);
T = length(epsi); 
Iepsi = zeros(length(epsi)) ;

if type == 1,
    for t = (maxi+1):T, 
        h(t) = para' * [1 ; epsi(t-(1:p)).^2;  h(t-(1:q))];
        epsi(t) = n(t)*sqrt(h(t));
        y(t) = [a ; b ; d]' * [1 ; y(t-(1:P)) ; epsi(t-(1:Q))] + epsi(t);
    end
elseif type == 2,
    if EST == 1,
        invtheta = 1/theta; 
        for t = (maxi+1):T, 
            h(t) = para' * [1 ; epsi(t-(1:p)).^2; expcdf(epsi(t-(1:p)),invtheta)*epsi(t-(1:p)).^2 ; h(t-(1:q))];
            epsi(t) = n(t)*sqrt(h(t));
            y(t) = [a ; b ; d]' * [1 ; y(t-(1:P)) ; epsi(t-(1:Q))] + epsi(t);
        end
    else
        F = zeros(T,1);
        F(1:maxi) = 1/(1+exp(-theta*epsi(1:maxi)));
        for t = (maxi+1):T,
            F(t) = 1/(1+exp(-theta*epsi(t-(1:p)))) - 0.5;
            h(t) = para' * [1 ; epsi(t-(1:p)).^2; F(t)*(epsi(t-(1:p)).^2);h(t-(1:q))];
            epsi(t) = n(t)*sqrt(h(t));
            y(t) = [a ; b ; d]' * [1 ; y(t-(1:P)) ; epsi(t-(1:Q))] + epsi(t);
        end
    end
elseif type == 3, 
    for t = (maxi+1):T,
        h(t) = para' * [1 ; epsi(t-(1:p)).^2;  h(t-(1:q)) ; Iepsi(t-(1:assym)).^2];
        epsi(t) = n(t)*sqrt(h(t));
        Iepsi(t) = (epsi(t)>0)*epsi(t); 
        y(t) = [a ; b ; d]' * [1 ; y(t-(1:P)) ; epsi(t-(1:Q))] + epsi(t);
    end
end

simulateEpsi = epsi((maxi+1+1000):T);
H = h(maxi+1+1000:T); 
Y = y(maxi+1+1000:T);

end
