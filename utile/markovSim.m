function states = markovSim(dim, transM)         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Simulation of a Markov Chain 
%
% INPUTS:
%   dim: length of the markov chain
%   transM: the transition matrix which have to be a square matrix
%
% OUTPUTS:
%   states: a matrix dim x k. The state of the Markov chain at any time. 
%               The rows represent the periods and the column the states.
%               For exemple, if the process is in state 1 at time t, the
%               cell (t,1) take value equal to 1 and cells (t,2:k) 0. 
%
% Author: Thomas Chuffart
% Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chekin' Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,
    error('You have forgotten some input');
end

if length(dim) > 1 || dim < 0 
    error('dim should be a positive scalar');
end

if isscalar(transM),
    error('transM has to be a squared matrix');
elseif size(transM,1) ~= size(transM,2),
    error('transM has to be a squared matrix');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% variables initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = size(transM,1);
T = dim;
states = zeros(T,k); 
states(1,1) = 1; 
prob = rand(T,1);

for i = 2:T,
    pastState=find(states(i-1,:)==1);
    if prob(i,1)<transM(pastState,pastState),
        states(i,pastState)=1;
    else    
        idxother=find(states(i-1,:)==0);
        prob2= transM(:,pastState);
        a=[transM(pastState,pastState) ; prob2(idxother)];
        cum_sum=cumsum(a);
        sorted=sort([cum_sum ; prob(i,1)]);
        idx=find(prob(i,1)==sorted)-1;      
        states(i,idxother(idx))=1;        
    end
end  