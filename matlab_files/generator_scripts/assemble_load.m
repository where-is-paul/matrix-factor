function [L]=assemble_load(p,t,te,edge_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembles the load vector for the mixed discretization of
% the lowest order. Boundary conditions are not taken
% into account.
% NOTE: Data has to be adjusted in data.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get point and element data 
%
el_num =length(t);                         % number of elements
x=p(1,:);                                  % global x-coordinates of vertices
y=p(2,:);                                  % global y-coordinates of vertices

L=zeros(edge_num,1);

% quadrature points on master-element
%
[xnod,ynod,weight]=quad5;

% number of quadrature points
%
quad_num = length(weight);

% get Nedelec reference shape functions[xnod,ynod,weight]=quad5;
%
sh_ref=Ned1_0(xnod,ynod);

for l=1:el_num                            % loop over elements
 
  % get global coordinates of vertices 1, 2, 3 of current element
  x0=x(t(1,l));   
  x1=x(t(2,l));
  x2=x(t(3,l));
  y0=y(t(1,l));
  y1=y(t(2,l));
  y2=y(t(3,l));
 
  J=[x1-x0,x2-x0;y1-y0,y2-y0];        % Jacobian of affine 
  trans=[x0;y0];                      % translation vector 
  JinvT=transpose(inv(J));            % transpose of inverse
  Jdet=det(J);                        % determinant of J
  
  for k=1:quad_num                    % loop over quadrature points
    
    dx=Jdet*weight(k);                % volume element  

    [quad]=J*[xnod(k);ynod(k)]+trans;       % physical point
    xquad=quad(1);  
    yquad=quad(2); 
    
    % get data at physical point
    f = data(xquad,yquad);
     
    for i=1:3    % loop over local dofs (index of test function v)
                       
      iglob_e=te(i,l);    % global dof-number of edge i of a
      
      %%%% Piola-transformation for i-th Nedelec shape-function
      %%%%% evaluated at quadrature point (x(k),y(k)) %%%%%
      
      shi=JinvT*sh_ref(:,i,k);      
   
      % assemble load-vector 
      % (pay attention to orientation te(i+3,l)=+-1 of edge!)
      
      L(iglob_e)=L(iglob_e) + te(i+3,l)*(f'*shi)*dx;  
      
    end
  end                                   
end

return








