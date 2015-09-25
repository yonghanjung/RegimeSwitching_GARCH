function [c ceq] = constrMSGARCH(x,k,nbcgtot)
%
% Purpose: function to compute the stationarity constraint of a MS-GARCH
%               model
%
% INPUTS: - x the vector of parameters 
%              - k the number of regimes
%              - nbcgtot: the number of GARCH parameters 
%
% OUTPUTS: - c inequality constraints
%
% Author: Thomas CHUFFART
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nbcg = nbcgtot/k;

    if k == 2,
        a1 = x(2);
        b1 = x(3);
        a2 = x(5);
        b2 = x(6);
        p = x(7);
        q = x(8);

        M = [p*(a1+b1) 0 (1-q)*(a1+b1) 0 ; p*a2 p*b2 (1-q)*a2 (1-q)*b2 ; (1-p)*b1 (1-p)*a1 q*b1 q*a1 ; 0 (1-p)*(a2+b2) 0 q*(a2+b2)];

    else
        beta = diag(x(3:nbcg:end-(k^2)));
        alpha = x(2:nbcg:end-(k^2));
        for i = 1:k,
            for j = 1:k,
                e = zeros(1,k);
                e(i) = 1; 
                A = x(nbcgtot+1+(i-1)+(k*(j-1))) *(beta + alpha*e);
                M(1+(k*(i-1)):1+(k*(i-1))+(k-1), 1+(k*(j-1)):1+(k*(j-1))+(k-1)) = A;
            end    
        end

    end
    [~,E] = eig(M);
    c = max(max(E))-1;
    ceq = [];   
end
    
    