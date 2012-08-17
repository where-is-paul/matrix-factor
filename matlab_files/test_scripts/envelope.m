function [ env ] = envelope( A )
%   calculates the envelope of a matrix
    env = 0;
    for i = 1:size(A,1)
        env = env + abs(find(A(i,:), 1, 'last')-i);
    end
end

