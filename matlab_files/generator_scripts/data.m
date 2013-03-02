%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provides data for the mixed discretization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function f = data_m(x,y)
%
% Anna S., 3.6.02
%
% gets RHS f at point (x,y).
%
% f(i): i-th component of RHS f at point(x,y).


f=zeros(2,1);

% % example 1   ( ==> div f = 0)
% %%%%%%%%%%%%%%%
 %f(1)=1;
 %f(2)=1;
% %%%%%%%%%%%%%%%

% example 2   ( ==> div f \neq 0)
%%%%%%%%%%%%%%%
f(1)=2; %-2.*x.*(1-y.^2);
f(2)=2; %-2.*y.*(1-x.^2);
%%%%%%%%%%%%%%%

% % example 3
% %%%%%%%%%%%%
% 
% f(1)=cos(pi.*x).*sin(pi.*y);
% f(2)=-sin(pi.*x).*cos(pi.*y);
% f= (2*pi*pi + c).*f;
% 
