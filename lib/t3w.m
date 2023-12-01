function [y] = t3w(M,x,transpose=0);
    % tensor3 by unwrapping the three operators
    % Take in transpose as a flag
    if transpose; 
        y = tensor3(M{3}',M{2}',M{1}',x);
    else;
        y = tensor3(M{3},M{2},M{1},x);  
    end;
end;