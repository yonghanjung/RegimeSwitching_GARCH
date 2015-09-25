function [thetahat results struct] = argarchEst(y,type,ORDERS,others,valdep,dist,hreal,options,printresults)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: 
%   This function estimates an univariate AR(P)-GARCH(p,q)
%   process under constraints.
%   You can choose the type of the GARCH process: GARCH, LST-GARCH or
%   GJR-GARCH.
%
% INPUTS:
%   y: The data of the process
%   type: The type of model to estimate. 1: GARCH 2: LST-GARCH and
%   3:GJR-GARCH
%   ORDERS: A vector of length = 5 with ORDERS(1) = C (=1 if constante in
%       the AR process), ORDERS(2) = P (AR order), ORDER(3) = Q (MA order), ORDER(4) = p (ARCH order)
%       and ORDER(5) = q (GARCH order)
%   others: If type 2, others define de smoothing function 
%              If type 3, others define the order of assymetry 
%   valdep: vector of initial values for the likelihood initialization 
%   dist: Distribution of the error term 'GAUSSIAN' or 'T' (Gaussian
%   available only for the moment)
%   options: A structure of options for fmincon 
%   printresults: A boolean, 1 if you want to print the results of the estimation
%                       when you call argarchEstim. 1 by defaults.
%
% OUTPUTS:
%   thetahat: Vector of estimate coefficient
%   results: a structure with many statistics (loss functions, information
%   criteria for example)
%   struct: a structure wich contain output of fmincon
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check the inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        c = ORDERS(1);
        P = ORDERS(2);
        Q = ORDERS(3);
        p = ORDERS(4);
        q = ORDERS(5);        
        
if nargin < 9, 
    printresults = 0;
end

if nargin < 8 || isempty(options), 
    options  =  optimset('fmincon');
    options  =  optimset(options , 'Algorithm ','sqp');
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
end

if nargin < 7,
    flagmse = 0;
else
    flagmse = 1;
    if size(hreal,2) > 1,
        error('hreal vector should be a column vector');
    end
end

if nargin < 3,
    if length(ORDERS) ~= 5,
        error('ORDERS has to be a vector of length 5');
    else
        c = ORDERS(1);
        P = ORDERS(2);
        Q = ORDERS(3);
        p = ORDERS(4);
        q = ORDERS(5);
        if (length(P) > 1) || P < 0 || length(Q) > 1|| Q < 0
            error('AR and MA should be positive scalars')
        end
        if (length(c) > 1) || c < 0 
            error('C should be positive scalar')
        end
        if (length(p) > 1) || (length(q) > 1) || p < 0 ||q < 0
            error('p and q should be positive scalars');
        end
        m = max(ORDERS);
    end
end

if nargin < 2,
    if isscalar(type),
        if type > 3 || type <1,
        error('type has to be equal at 1, 2 or 3');
        end
    else
        error('type has to be a scalar');
    end
end



if nargin < 6,
    fprintf('Dist is GAUSSIAN by default');
elseif strcmp(dist,'GAUSSIAN')  == 0 && strcmp(dist,'T') == 0,
    error('bad distribution');
end 

if nargin < 5 || isempty(valdep),
    valdep = 0;
end

if nargin < 4 || isempty(others),
    fprintf('\n Other is empty or not call. if GARCH, no problem, if ST-GARCH, others is LOGIT by default, if GJR others = 1.');
    if type == 1,
        others = 0;
        assym = 0;
    elseif type ==2,
        assym = 0; 
        others = 'L';
    elseif type == 3,
        assym = 1;
        others = assym;
    end
else
    if type == 2, 
        if ischar(others) ,
            EST = strcmp('E',others) ;
            LST = strcmp('L',others) ;
            if EST == 0  &&  LST == 0, 
                error('Erreur in smooth function input, has to be equal to a string variable L or E');
            end
        else
            others = 'L'; 
        end
    elseif type == 3,
        if isscalar(others),
            assym = others; 
        else
            error('Error, others is not a scalar');
        end
    end
end

if p == 0,
    q = 0; 
    nbcoef = P + c + Q; 
    nbcoefgarch = 0;
    nbcoefar = nbcoef;
elseif type == 1, 
    nbcoef = sum(ORDERS) + 1; 
    nbcoefgarch = 1 + p + q ; 
    nbcoefar = Q + P + c;
elseif type == 2, 
    nbcoef = sum(ORDERS) + 2 + p; 
    nbcoefgarch = 1 + (2*p) + q; 
    nbcoefar = Q + P + c;
elseif type == 3, 
    nbcoef =  sum(ORDERS)+ 1 + assym;
    nbcoefgarch = 1 + p + q + assym; 
    nbcoefar = Q + P + c;
end

if length(valdep) ~= nbcoef,
    if valdep == 0,
    fprintf('\n Valdep is empty, we assign it default value \n'); 
    else
        fprintf('\n valdep is mispecified, we assign it default values \n')
    end
    
    if P ~= 0 || c ~= 0,
            [Y X] = matrixlag(y,P,c);
            bhat = X\Y; 
            bhat = real(bhat);
    else
        bhat = [] ;
    end
    
    if Q ~= 0,
        bhat = [bhat ; ones(Q,1)/2];
    end
    
    if type == 1,   
        alpha  =  .05*ones(p,1);
        beta   =  0.05*ones(q,1)/q;
        omega  =  0.01;  
        if isempty(bhat),
            valdep = [omega ; alpha ; beta];
        else
            if p ~= 0,
                valdep = [bhat ; omega ; alpha ; beta];
            else
                valdep = bhat;
            end
        end
    elseif type == 2,
        alpha  =  .05*ones(2*p,1);
        beta   =  0.05*ones(q,1)/q;
        omega  =  0.01;
        theta = 2;
        if isempty(bhat),
            valdep = [omega ; alpha ; beta ; theta];
        else
            if p ~= 0,
                valdep = [bhat ; omega ; alpha ; beta ; theta];
            else
                valdep = bhat;
            end
        end
        
    elseif type == 3
        alpha  =  .05*ones(p,1);
        coeffassym  =  0.05*ones(assym,1)/assym;
        beta   =  0.05*ones(q,1)/q;
        omega  =  0.01;
        if isempty(bhat),
            valdep = [omega ; alpha ; beta ; coeffassym];
        else         
            if p ~= 0,
                valdep = [bhat ; omega ; alpha ; beta ; coeffassym];
            else
                valdep = bhat;
            end
        end
    end
end

if size(valdep,2)>1
    valdep = valdep';
end

if nargin < 3,
    error('Wrong input argument');
end

if size(y,2) > 1 
   error('Data vector should be a column vector');
end

data = y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate the parameters and compute the standard erreur:
% Constraint on the parameters:

if type == 1
    A = [zeros(nbcoefgarch,nbcoefar) -eye(nbcoefgarch) ; zeros(1,nbcoefar+1) ones(1,nbcoefgarch-1) ];
    bt = [zeros(1, nbcoefgarch) 1 - (1e-6)];
    LB  = [-ones(nbcoefar,1) ; zeros(nbcoefgarch,1)];
    UB  = [];
elseif type ==2,
    nbcontr = (2*(P+Q))+(3*p)+q+2;
    para = nbcoef; 
    A = zeros(nbcontr,para);
    flag = 1;
        k = 1;
        while flag == 1,
            if P == 0,
                flag = 2;
                break
            end
            for j = 1:P+Q,
                A(k,j+c) = -1;
            end
            if k == (2*(P+Q)),
                flag = 2;
            end
            k = k + 1;
        end
        A(k,(2*(P+Q))+c+1) = -1;
        k = k + 1;
        j = 1;
        while flag == 2,
            A(k,P+Q+c+1+j) = -1;
            j = j + 1;
            if k == 2*(P+Q)+1+p,
                flag = 3;
            end
            k = k + 1;
        end
        j=1;
        while flag == 3,
            A(k,P+Q+c+1+j) = -2;
            A(k,P+Q+c+1+j+p) = -1;
            if k == (2*(P+Q))+1+(2*p),
                flag = 4;
            end
            k = k + 1;
        end
        j = 1;
       while flag == 4,
           A(k,P+Q+c+1+j) = -2;
           A(k,P+Q+c+1+j+p) = 1;
           if k == (2*(P+Q))+1+(3*p),
               flag = 5;
           end
           k = k + 1;
       end
       j = 1;
       while flag == 5,
           A(k,P+Q+c+1+(2*p)+j) = -1;
           if k ==  (2*(P+Q))+1+(3*p)+q, 
               flag = 6;
           end
           j = j +1;
           k = k+1;
       end        

    A(end, end) = -1;  
    theta = valdep(end);
    bt = zeros(1,nbcontr); 
    LB = [0, -ones(1,nbcoef-2) 0.0001];
    UB=[1, ones(1,nbcoef-2) (theta+10*theta)];
elseif type == 3,
    A = [zeros(nbcoefgarch-1,nbcoefar) -eye(nbcoefgarch-1) ; zeros(1,nbcoefar+1) ones(1,nbcoefgarch-2) ];
    m = size(A); 
    A = [A zeros(m,1)];
    bt = [zeros(1, nbcoefgarch-1) 1 - (1e-6)];
    LB  = [-ones(nbcoefar,1) ; zeros(nbcoefgarch-1,1) ; -ones(assym,1)];
    UB  = [];
end

% We call fmincon function for the estimation

[thetahat,fval,exitflag,output,lambda,grad,HESSIAN]  = fmincon('argarchLik',valdep,A,bt,[],[], ...
                                                            LB,UB,[],options,y,ORDERS,others,nbcoef,type); 
[sumLLF,LLF,ht,epsi] = argarchLik(thetahat,y,ORDERS,others,nbcoef,type);
struct.output = output;
struct.fval = fval;
struct.exit = exitflag;
struct.lamb = lambda;
struct.grad = grad;
struct.hessian= HESSIAN;

%Compute robust standard errors with the sandwich matrix

[VCV]=robustvcv(@argarchLik,thetahat,0,y,ORDERS,others,nbcoef,type);
stderr=sqrt(diag(VCV));

T = length(epsi);
if P == 0 && c == 0,
    y = epsi;
elseif P == 0,
    y = thetahat(1) + epsi; 
else
    if c == 0,
        thetahat = [0 ; thetahat];
    end
    for t = m+1:T,
        y(t) = thetahat(1:c+P+Q)'*[1:c;y(t-(1:P));epsi(t-(1:Q))'] + epsi(t);
    end
end

%We group all our results in a structure

results.sumLLF = sumLLF;
results.stderr = stderr;
results.LLF = LLF;
results.AIC = InfoCrit(-sumLLF,thetahat,y,1);
results.BIC = InfoCrit(-sumLLF,thetahat,y,2);
results.H = ht;
results.epsi = epsi;

if flagmse,
    dim = length(y);
    results.mise = lossfun(hreal(1:dim),ht,2);
    results.qlike = lossfun(hreal(1:dim),ht,3);
    results.mae = lossfun(hreal(1:dim),ht,4);
    
    results.mise1 = lossfun(data(1:dim).^2,ht,2);
    results.qlike1 = lossfun(data(1:dim).^2,ht,3);
    results.mae1 = lossfun(data(1:dim).^2,ht,4);
end


if type == 0,
    % Affichage
    if printresults,
        B = strcat(char('C',...
        char(regexprep(regexp(sprintf('AR%1.0f/',0:P), '/', 'split'), '^.*0', '')),...
        'K', ...
        char(regexprep(regexp(sprintf('ARCH%1.0f/',0:p), '/', 'split'), '^.*0', '')),...
        char(regexprep(regexp(sprintf('GARCH%1.0f/',0:q), '/', 'split'), '^.*0', '')),...
        char(regexprep(regexp(sprintf('GJR%1.0f/',0:assym), '/', 'split'), '^.*0', ''))));

        fprintf('-------------------------------------------------\n');
        fprintf('Convergence achieved after %1.0f iterations\n', OUTPUT.iterations);
        fprintf('-------------------------------------------------\n');
        fprintf('Parameters  Coefficients  Std Errors    T-stats   p-val\n');
        fprintf('-------------------------------------------------\n');
        for i = 1:size(thetahat,1)
            if thetahat(i) < 0
                fprintf(strcat('  %s     %1.',num2str(round(5-length(sprintf('%1.0f', abs(thetahat(i)))))),...
                'f %1.',num2str(round(5-length(sprintf('%1.0f', stderr(i))))),'f   %1.',...
                num2str(round(5-length(sprintf('%1.0f',abs(resultsSTTest.tstat(i)))))),'f\n'),...
                'f %1.',num2str(round(5-length(sprintf('%1.0f', resultsSTTest.pval(i))))), B(i,:),...
                thetahat(i), stderr(i), resultsSTTest.tstat(i), resultsSTTest.pval(i));
            else    
                fprintf(strcat('  %s      %1.',num2str(round(5-length(sprintf('%1.0f', abs(thetahat(i)))))),...
                'f %1.',num2str(round(5-length(sprintf('%1.0f', stderr(i))))),'f        %1.',...
                num2str(round(5-length(sprintf('%1.0f', abs(resultsSTTest.tstat(i)))))), 'f %1.',...
                num2str(round(5-length(sprintf('%1.0f', resultsSTTest.pval(i))))),'f\n') , ...
                B(i,:), thetahat(i), stderr(i), resultsSTTest.tstat(i), resultsSTTest.pval(i));
            end
        end
        fprintf('-------------------------------------------------\n');
        fprintf('R-Squared: %1.4f\n', results.Rsquared);
        fprintf('Adjusted R-Squared: %1.4f\n', results.AdjRsquared);
        fprintf('Log Likelihood: %1.0f\n', LLF);
        fprintf('Akaike Information Criteron: %1.0f\n', results.AIC);
        fprintf('Bayesian Information Criteron: %1.0f\n', results.BIC);
        fprintf('Kurtosis des résidus standardiséss: %1.4f\n', results.kurt);
        fprintf('Skewness des résidus standaridisés: %1.4f\n', results.skew);
        fprintf('Jarque-Bera statistic and p-value: %1.4f %1.4f %1.0f \n\n', results.JBstat, results.JBpval);
        fprintf('-------------------------------------------------\n');

    end
end
end

