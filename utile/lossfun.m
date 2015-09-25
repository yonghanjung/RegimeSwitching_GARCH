function f = lossfun(real,predict,flag,b)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PURPOSE: this function compute several loss functions like MSE, MISE,
%                 quadratic loss function...
%
% INPUTS:
%   - real: real or value of the paramater you study, for exemple ht can be
%             the realized volatility 
%   - predict: the predict value of the parameter of interest hat(ht) for
%               exemple, with ht wich is a type garch process
%   - flag: 1 for mean squared erro
%             2 for mean integrated squared error
%
% OUPUT:
%   - f, the value of tge loss function. 
%
% Author: Thomas CHUFFART
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3,
    error('wrong number of arguments');
end

if size(real,1) ~= size(predict,1) || size(real,2) ~= size(predict,2),
    error('real and predict have to be of same size');
end

switch flag,
    case 1 ,
        f = (real - predict).^2;
    case 2
        dim = length(predict);
        f = sum((real-predict).^2)/dim;
    case 3
        dim = length(predict);
        f = sum((log(real)+predict./real))/dim; 
    case 4
        dim = length(predict);
        f = sum(abs(real - predict))/dim;
end