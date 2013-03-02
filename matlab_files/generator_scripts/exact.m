function [u, curlu, p, gradp] = exact(x,y)
%
% Anna Schneebeli, June 02
%
% 
% u(i,k): i-th component of exact solution u at point x(k), y(k)
% curlu(k) : curl of exact solution at point x(k), y(k)
% p(k): exakt p at point x(k), y(k)
% gradp(i,k): i-th component of gradient of p at point x(k), y(k)

u=zeros(2,length(x));

% % example 1  
% %%%%%%%%%%%%%%%%
% u(1,:)=1-y.^2;
% u(2,:)=1-x.^2;
% 
% curlu = 2.*y-2.*x;
% 
% p = 0;
% 
% gradp(1,:)=0;
% gradp(2,:)=0;
% 
% %%%%%%%%%%%%%%%%%

% example 2
%%%%%%%%%%%%%%%%
u(1,:)=1-y.^2;
u(2,:)=1-x.^2;

curlu = 2.*y-2.*x;

p = (1-x.^2)*(1-y.^2);

gradp(1,:)=-2.*x.*(1-y.^2);
gradp(2,:)=-2.*y.*(1-x.^2);

%%%%%%%%%%%%%%%%%

% % example 3
% %%%%%%%%%%%%%%%%%
% u(1,:)=cos(pi.*x).*sin(pi.*y);
% u(2,:)=-sin(pi.*x).*cos(pi.*y);
% 
% curlu = -2*pi.*cos(pi.*x).*cos(pi.*y);
% %%%%%%%%%%%%%%%%%%






