function [rho, acfi , stdfi] = autocorr(X,k,graph)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose:
%   This function computes the empirical auto-correlation.
%
% INPUTS:
%   - X: The data
%   - k: The number of lag
%
% OUTPUTS:
%   - rho: the empirical autocorrelation 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Checkin' INPUTS%%%%%%%%%%%%%%%%%%%
if nargin < 2,
    error('2 inputs are needed');
end

if nargin<3,
    graph = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = zeros(k,1); 
T = length(X);
for i = 1:k,
    rho(i) = autocov(X,i)/std(X)/std(X(i+1:T));
end

%compute the robust standard deviation: 

if graph == 1,

    for i = 1:k,
        Y = X(i+1:T);
        [~, Xlag] = matrixlag(X,1,1);
        Xlag = Xlag(1:T-i,:);
        dim=length(Y);
        A = ((Xlag'*Xlag)/length(Y))^(-1);
        v = repmat(Y-mean(Y),1,2); 
        S = ((Xlag.*v)'*(Xlag.*v))/length(Y);
        vc = (A*S*A)/length(Y);
        stdrb(i) = sqrt(vc(2,2));
        stderr(i)=sqrt(1/dim);
    end
    
    stdrb = stdrb';
    stderr = stderr';
    
    acfi = bar(rho); hold on;
    stdfi = plot((0:k+1)', [1.96 *[stdrb(1); stdrb ; stdrb(k)] -1.96*[stdrb(1); stdrb ; stdrb(k)]], '-.','Color','r'); 
    stdfi2 = plot((0:k+1)', [1.96 *[stderr(1); stderr ; stderr(k)] -1.96*[stderr(1); stderr ; stderr(k)]],'-.','Color','b');
    axis auto; 
    t = title('Autocorrelation and Robust standard error');
    set(t,'FontSize',10)
    colormap summer;
    hold off;
end

end