function mat=m2d(beta,gamma,n);
h=1./(n+1); 
a=4; b=-1-gamma; c=-1-beta; d=-1+beta; e=-1+gamma; 
I=speye(n);
t1=spdiags([c*ones(n,1),2*ones(n,1),d*ones(n,1)],-1:1,n,n);
t2=spdiags([b*ones(n,1),2*ones(n,1),e*ones(n,1)],-1:1,n,n);
mat=kron(I,t1)+kron(t2,I);