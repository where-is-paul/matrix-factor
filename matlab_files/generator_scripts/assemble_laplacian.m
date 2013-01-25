function [S]=assemble_scalar_lapacian(p,t)
%
% Assemble the mass and stiffness matrix of the Laplacian 
% with piecewise linear shape functions
% No boundary conditions are applied

disp('assembling Laplacian....')

% get point and element data 

el_num =length(t);                         % number of elements
vertex_num = length(p);                    % number of vertices
x=p(1,:);                                  % physical x-coordinates of vertices
y=p(2,:);                                  % physical y-coordinates of vertices

S=sparse(vertex_num,vertex_num);

% gradients of nodal reference shape functions:

ref_nodal_grads = [[-1;-1],[1;0],[0;1]];


for l=1:el_num          % loop over all elements                        

  % get global coordinates of vertices 1, 2, 3 of current element
  x0=x(t(1,l));   
  x1=x(t(2,l));
  x2=x(t(3,l));
  y0=y(t(1,l));
  y1=y(t(2,l));
  y2=y(t(3,l));
 
  J=[x1-x0,x2-x0;y1-y0,y2-y0];          % Jacobian of affine element-map      
  trans=[x0;y0];                        % translation vector of affine element-map
  JinvT=transpose(inv(J));              % transpose of inverse
  Jdet=det(J);                          % determinant of J
  area_elem=Jdet*0.5;                   % area of the element  
  nodal_grads = JinvT*ref_nodal_grads;  % transformed gradients of nodal shape functions
  nodal_grads = JinvT*ref_nodal_grads;
       
  for i=1:3            % loop over local dofs (index of test function v)
    iglob_v=t(i,l);    % global dof-number of vertex i of actual element

    for j=1:3          % loop over local dofs (index of ansatz function v)     
      jglob_v=t(j,l);  % global dof-number of vertex j of actual element 

       % assemble stiffness matrix S

       dS=transpose(nodal_grads(:,i))*nodal_grads(:,j)*area_elem;
       S(iglob_v,jglob_v)=S(iglob_v,jglob_v)+dS;

       % assemble mass matrix M



       M(iglob_v,jglob_v)=0;
       
    end
  end                                   
end


return



