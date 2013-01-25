function [dxsh,dysh]=dshap1(x,y)
%
% computes derivatives of P1 shape functions on
% reference triangle on a list of points
% 
% dxsh(i,j): d/dx of shape function i on point x(j), y(j)
% dysh(i,j): d/dy of shape function i on point x(j), y(j)

dxsh=zeros(3,length(x));
dysh=zeros(3,length(y));
dxsh(1,:)=-1;
dxsh(2,:)=1;
dxsh(3,:)=0;
dysh(1,:)=-1;
dysh(2,:)=0;
dysh(3,:)=1;


