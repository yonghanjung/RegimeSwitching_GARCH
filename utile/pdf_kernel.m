function f = pdf_kernel(arg, sample, h )
%
% Purpose: compute the non parametric estimation of a density function by
%               kernel method
%
%  INPUT: - arg: the argument we want to estimate the density function
%             - sample: the support of the estimate function
%             - h: the brandwith parameter
%
%  OUTPUT: - f the estimation of the density funcion
%
%  Author: Thomas CHUFFART
%  Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check'in input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = length (arg);
f = zeros (dim,1);
for i=1:dim,
    loc = (arg(i)-sample)/h;
    f(i) = mean(dnorm(loc, 0,1));
end
f = f/h;
end

