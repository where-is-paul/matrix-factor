%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
k=0;
gamma=1;
for i=1:2,
    load(['matrices',num2str(i),'.mat']); 
    [m,n]=size(B);
    F=A-k^2*M;
    Q=[F B'; B sparse(m,m)];
    P =[A+gamma*M  sparse(n,m); sparse(m,n) L];
    PL=[F+B'*(L\B) sparse(n,m); sparse(m,n) L];
    PW=[F+B'*(W\B) sparse(n,m); sparse(m,n) W];
    PI=[F+B'*B sparse(n,m); sparse(m,n) speye(m)];
    eig_P  = sort(real(eig(full(Q),full(P))));
    eig_PL = sort(real(eig(full(Q),full(PL))));
    eig_PW = sort(real(eig(full(Q),full(PW))));
    eig_PI = sort(real(eig(full(Q),full(PI))));
    save(['preconditioners',num2str(i),'.mat']); 
end








