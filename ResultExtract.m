function [Result] = ResultExtract(Smooth)
    Result = [];
    [n,m] = size(Smooth);
    for idx = 1 : n,
        if Smooth(idx,1) > Smooth(idx,2),
            Result = [Result;0];
        else
            Result = [Result;1];
        end
    end
    
          
            