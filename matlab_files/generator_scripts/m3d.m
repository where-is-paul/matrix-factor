function mat=m3d(beta,gamma,delta,n);
h=1./(n+1); 
a=6; b=-1-gamma; c=-1-beta; d=-1+beta; e=-1+gamma; f=-1-delta; g=-1+delta; 
I=speye(n);
t1=spdiags([c*ones(n,1),a*ones(n,1),d*ones(n,1)],-1:1,n,n);
t2=spdiags([b*ones(n,1),zeros(n,1),e*ones(n,1)],-1:1,n,n);
t3=spdiags([f*ones(n,1),zeros(n,1),g*ones(n,1)],-1:1,n,n);
%mat=kron(I,t1)+kron(t2,I);
mat=kron(I,kron(I,t1))+kron(I,kron(t2,I))+kron(t3,kron(I,I));