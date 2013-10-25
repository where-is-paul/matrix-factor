function [L, D, p, S, B] = ildl(A, opts)
% version of ildl that applies a shift matrix to the specified
% problem. to be more precise, the block diagonal matrix  with
% diagonal block alpha * ||A||_1 * [0 1; 1 0] is added to the 
% original matrix, producing a positive shift in its eigenvalues.
    
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
        opts.pp_tol = 2.00;
    end
    if (~any(strcmp('ordering',fieldnames(opts))))
        opts.ordering = 'amd';
    end
    if (~any(strcmp('shift_factor',fieldnames(opts))))
        opts.shift_factor = 0.000;
    end
    if (~any(strcmp('equil',fieldnames(opts))))
        opts.equil = 'y';
    end
    
    % uncomment to display given options.
    display(opts);
    
    if (~issparse(A))
        error('Input matrix must be sparse');
    end
    if (norm(A-A',1) ~= 0)
        error('Input matrix must be symmetric');
    end
    % sort of wasteful, but we'll perform a first factorization
    % to get A equilibriated and reordered (also with BKP off)
    % and then a second factorization on the B, with equilibriation
    % and permutation turned off 
    [~, ~, p1, S, C] = ildl(A, opts.fill_factor, opts.tol, ...
                            0.00, opts.ordering, opts.equil);  
    
    % turn off equilibriation and ordering
    opts.equil = 'n';
    opts.ordering = 'none';
    
    display('Condition number before shifting');
    display(condest(C));
    % keyboard;
    % perform a shift on our matrix
    C = C + block_shift(size(A,1), opts.shift_factor);
    
    display('Condition number after shifting');
    display(condest(C));
    %{
    n = size(A,1);
    d = mod(1:n,2);
    shift = spdiags(d', -1, n, n) + spdiags(d', -1, n, n)';
    if (mod(n,2) == 1)
        shift(n,n) = 1;
    end
    C = C+opts.shift_factor*shift;
    %}

    [L, D, p2, ~, B] = ildl(C, opts.fill_factor, opts.tol, ...
                       opts.pp_tol, opts.ordering, opts.equil);
    p = p1(p2, :);
        
end

