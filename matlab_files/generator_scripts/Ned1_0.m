function shapef=Ned1_0(x,y)
%
% Anna Schneebeli, 3.6.02
% 
% computes first type Nedelec-shape functions of lowest dergee on
% reference triangle on a list of points
% 
% shapef(i,j,k): i-th component of shape function j at point x(k), y(k)
%

shapef=zeros(2,3,length(x));

shapef(1,1,:)=1-y;
shapef(2,1,:)=x;

shapef(1,2,:)=-y;
shapef(2,2,:)=x;

shapef(1,3,:)=-y;
shapef(2,3,:)=x-1;

return