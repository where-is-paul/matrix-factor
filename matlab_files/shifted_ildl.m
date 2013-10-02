function [L, D, p, S, B] = ildl(A, opts)
% version of ildl that applies a shift matrix to the specified
% problem. to be more precise, the block diagonal matrix  with
% diagonal block alpha * [0 1; 1 0] is added to the original matrix,
% producing a positive shift in its eigenvalues.
    
    if (nargin < 2)
        opts.type = 'default';
    end
    % default fill factor
    if (~any(strcmp('fill_factor',fieldnames(opts))))
        opts.fill_factor = 2.0;
    end
    if (~any(strcmp('tol',fieldnames(opts))))
        opts.tol = 0.001;
    end
    if (~any(strcmp('pp_tol',fieldnames(opts))))
        opts.pp_tol = 1.00;
    end
    if (~any(strcmp('ordering',fieldnames(opts))))
        opts.ordering = 'amd';
    end
    if (~any(strcmp('shift_factor',fieldnames(opts))))
        opts.shift_factor = 0.001;
    end
    
    display(opts);
    
    n = size(A,1);
    d = mod(1:n,2);
    shift = spdiags(d', -1, n, n) + spdiags(d', -1, n, n)';
    A = opts.shift_factor*shift;
    [L D p S B] = ildl(A, opts.fill_factor, ...
                       opts.tol, opts.pp_tol, opts.ordering);
    
    
end

