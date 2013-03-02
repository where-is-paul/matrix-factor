%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=0

for i=1:7
    i
    disp('Loading matrices...');
    load(['matrices_homogeneous',num2str(i),'.mat']); 
    
    [m,n]=size(B);
    Q=[A-kk^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];

    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    
    disp('Running PCG');
    [xcg,flagcg,relrescg,itercg,resveccg] = mypcg(A-kk^2*M,f,1e-10,1000,A+(1-kk^2)*M,[],[],L,B); flagcg,itercg
    
    
    %AW=A+B'*(W\B);
    
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AL,f,1e-12,100,A+M); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AW,f,1e-12,100,A+nel*B'*B); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AI,f,1e-10,100,A+); ITER,FLAG
    %pause
    %disp('Setting up preconditioners...');
    P =[A+(1-kk^2)*M  sparse(n,m); sparse(m,n) L];
    %P2 =[A+M  sparse(n,m); sparse(m,n) W];
    %PL=[AL-kk^2*M sparse(n,m); sparse(m,n) L];
    %PW=[AW-kk^2*M sparse(n,m); sparse(m,n) W];
    %I=speye(m);
    %delta=1/nel;
    %PI=[A-kk^2*M+1/delta*B'*B sparse(n,m);sparse(m,n) delta*I];
    
    %disp('Computing eigenvalues...');
    %eig_P  = sort(real(eig(full(Q),full(P)))),pause
   %eig_PL = sort(real(eig(full(Q),full(PL)))),pause
    %eig_PW = sort(real(eig(full(Q),full(PW)))),pause
    %eig_P2 = sort(real(eig(full(Q),full(P2)))),pause
    %eig_PI = sort(real(eig(full(Q),full(PI))));
    %disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); flag, iter
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,PI); flag, iter
    %[x2,flag2,relres2,iter2,resvec2,resveccg2] = minres(Q,rhs,1e-11,500,P2); flag2, iter2
    
    %[xL,flagL,relresL,iterL,resvecL,resveccgL] = minres(Q,rhs,1e-10,500,PL); flagL, iterL
    %[xW,flagW,relresW,iterW,resvecW,resveccgW] = minres(Q,rhs,1e-10,500,PW); flagW, iterW
    
    %pause
    %DUMMY = sparse(m+n,m+n);
    %tic,[xW2,flagW2,relresW2,iterW2,resvecW2,resveccgW2] = myminres(Q,rhs,1e-10,500,DUMMY,[],[],L,A,B); flagW2,iterW2,toc 
    
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








