function B = blokdiag(varargin)
%
% Purpose:
%   use: out = blockdiag(a,m)
%   This function returns a matrix B = [a 0 0 ; 0 a 0 ; 0 0 a] with a a
%   numeric elements (scalar, vector or matrix), m is the number of a on
%   the diagonal. In the exemple, m = 3. The zeros have the same size of a,
%   for exemple if a = [1 1 ; 1 1] an m = 2, B = [1 1 0 0]
%                                                                  |1 1 0 0|
%                                                                  |0 0 1 1|
%                                                                  |0 0 1 1|
%                                                                 
%   INPUT: varaging, with A a matrix and m a scalar
%
%   OUTPUT: B, the block diagonal matrix of A
%
%  Author: Thomas CHUFFART
%  Mail: thomas.chuffart@univ-amu.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check the inputs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isscalar(varargin{2}) == 0
        error('The second argument has to be a scalar');
    end

% Program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    varargin = repmat(varargin(1),[1,varargin{2}]);

    m = length(varargin);

    if m == 1, B = varargin{1}; return, end      
    Z = zeros(1,length(varargin{1}));% Get the sizes of each block

    for i = 1:m  
      s(i,:) = size(varargin{i});
    end
    % Build B one row at a time.
    for i = 1:m  
      Zl = repmat(Z,s(i,1),i-1);     % left zeros
      Zr = repmat(Z,s(i,1),m-i);  % right zeros
      row = [ Zl varargin{i} Zr];                % build the row
      if i ==1,
          B = row;
      else
          B = [B; row];         
      end
    end

end