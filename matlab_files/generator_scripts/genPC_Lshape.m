%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=1/8;

for ii=1:5
    ii
    disp('Loading matrices...');
    load(['Lshape_matrices',num2str(ii),'.mat']); 
    [m,n]=size(B);
    Q=[A-kk^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];

    
    disp('Running PCG');
    [xcg,flagcg,relrescg,itercg,resveccg] = mypcg(A-kk^2*M,f,1e-10,1000,A+(1-kk^2)*M,[],[],L,B); flagcg,itercg
    
    
    
    %condest(L),pause
    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    
    %condest(AL)
    %condest(AW),pause
    
    %AW=A+B'*(W\B);
    %W=normest(A)/normest(B'*B)*speye(size(B,1));
    %nAB=norm(A,1)/(norm(B,1)^2);
    %AW=A+B'*(W\B);
    
    
    
    %disp('Setting up preconditioners...');
    %P =[A+(1-kk^2)*M  sparse(n,m); sparse(m,n) L];
    %PL=[AL-k^2*M sparse(n,m); sparse(m,n) L];
    %PW=[AW-k^2*M sparse(n,m); sparse(m,n) W];
    
    %I=speye(m);
    %delta=1/nel;
    %PI=[A-kk^2*M+B'*B sparse(n,m);sparse(m,n) I];
    
    %[xcg,flagcg,relrescg,itercg,resveccg] = PCG(AL,f,1e-10,1000,A+(1-kk^2)*M); flagcg,itercg
    
    %condest(PL),pause
    
    %sort(eig(full(Q),full(PW))),pause
    
    %disp('Computing eigenvalues...');
    %eig_P  = sort(real(eig(full(Q),full(P))));
    %eig_PL = sort(real(eig(full(Q),full(PL))));
    %eig_PW = sort(real(eig(full(Q),full(PW))));
    %eig_PI = sort(real(eig(full(Q),full(PI))));
    %disp('Solving with MINRES...')
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); flag, iter
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,PI); flag, iter
    %[xL,flagL,relresL,iterL,resvecL,resveccgL] = minres(Q,rhs,1e-6,500,PL); flagL, iterL
    %tic,[xW,flagW,relresW,iterW,resvecW,resveccgW] = minres(Q,rhs,1e-10,500,PW); flagW, iterW,toc
    
    %DUMMY = sparse(m+n,m+n);
    %tic,[xL2,flagL2,relresL2,iterL2,resvecL2,resveccgL2] = myminres(Q,rhs,1e-8,500,DUMMY,[],[],L,A,B); flagL2,iterL2,toc 
    
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








