function [sh]=shap1(x,y)
%
% computes P1 shape functions on
% reference triangle on a list of points
% 
% sh(i,j): shape function i on point x(j), y(j)
%

sh=zeros(3,length(x));
sh(1,:)=1-x-y;
sh(2,:)=x;
sh(3,:)=y;


