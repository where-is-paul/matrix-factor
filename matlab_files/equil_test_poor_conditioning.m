opts.fill_factor = 2.5;
opts.tol = 0.0001;
opts.pp_tol = 1.0;
opts.ordering = 'amd';


iters = 10;
for k = 1:iters
    A = sprandsym(10000, 0.0005, 1e-10);
    A = A + 1e-8*speye(size(A,1));
    
    opts.equil = 'bunch';
    [l, d, p, S, B] = ildl(A, opts.fill_factor, opts.tol, ...
                       opts.pp_tol, opts.ordering, opts.equil);
   
    fprintf('The relative residual is %f.\n', norm(B - l*d*l', 1)/norm(B, 1));
    fprintf('The fill factor is %.3f.\n', nnz(l+d+l')/nnz(B));
    fprintf('The largest elem. of L is %.3f.\n', full(max(max(abs(l)))));
	fprintf('A has %i nnz.\n', nnz(A));
    fprintf('The condition number is %e.\n', condest(B));
    
    e = ones(size(B,1),1);
    b = S^(-1)*e;
    %[y, flag, relres, iter, resvec] = ...
        ld = @(x) d\(l\x);
        gmres(B,b,min(100,size(B,1)),1e-8,10,ld,l');

    fprintf('============================================\n');
    
    opts.equil = 'iter';
    [l, d, p, S, B] = ildl(A, opts.fill_factor, opts.tol, ...
                       opts.pp_tol, opts.ordering, opts.equil);
   
    fprintf('The relative residual is %f.\n', norm(B - l*d*l', 1)/norm(B, 1));
    fprintf('The fill factor is %.3f.\n', nnz(l+d+l')/nnz(B));
    fprintf('The largest elem. of L is %.3f.\n', full(max(max(abs(l)))));
	fprintf('A has %i nnz.\n', nnz(A));
    fprintf('The condition number is %e.\n', condest(B));
    
    e = ones(size(B,1),1);
    b = S^(-1)*e;
    %[y, flag, relres, iter, resvec] = ...
        ld = @(x) d\(l\x);
        gmres(B,b,min(100,size(B,1)),1e-8,10,ld,l');
        
    fprintf('\n');

end