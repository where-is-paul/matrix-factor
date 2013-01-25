%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1/8;
for i=3:3,
    i
    disp('Loading matrices...');
    load(['matrices',num2str(i),'.mat']); 
    [m,n]=size(B);
    disp('[m; n; n+m]');
    [m; n; m+n]
    Q=[A-k^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];
    disp('Setting up preconditioners...');
    P =[A+M  sparse(n,m); sparse(m,n) L];

    disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); 
    flag, iter,relres,resvec
 
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








