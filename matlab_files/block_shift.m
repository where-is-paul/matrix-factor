function [ C ] = block_shift( n, shift_factor )
% produces a shift matrix of size nxn
% we shift each 2x2 block on the diagonal of A
% by [0 1; 1 0]. if A has odd dimension, we shift
% A(n,n) by 1 to compensate for the parity.

d = mod(1:n,2);
shift = spdiags(d', -1, n, n) + spdiags(d', -1, n, n)';
if (mod(n,2) == 1)
    shift(n,n) = 1;
end
C = shift_factor*shift;


end

