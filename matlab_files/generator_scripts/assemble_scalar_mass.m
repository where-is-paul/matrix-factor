function [M]=assemble_scalar_mass(p,t)

%
% Assemble the mass matrix of the Laplacian 
% with piecewise linear shape functions
% No boundary conditions are applied

% get point and element data 

el_num =length(t);                         % number of elements
vertex_num = length(p);                    % number of vertices
x=p(1,:);                                  % physical x-coordinates of vertices
y=p(2,:);                                  % physical y-coordinates of vertices

M=sparse(vertex_num,vertex_num);

% quadrature points on master-element
%
[xnod,ynod,weight]=quad5;
quad_num = length(weight);

% nodal shape functions on quadrature points
%
sh_ref=shap1(xnod,ynod);



for l=1:el_num          % loop over all elements                        

  % get global coordinates of vertices 1, 2, 3 of current element
  x0=x(t(1,l));   
  x1=x(t(2,l));
  x2=x(t(3,l));
  y0=y(t(1,l));
  y1=y(t(2,l));
  y2=y(t(3,l));
 
  J=[x1-x0,x2-x0;y1-y0,y2-y0];          % Jacobian of affine element-map      
  Jdet=det(J);                          % determinant of J
  
  for k=1:quad_num  
    for i=1:3             % loop over local dofs (index of test function v)
      iglob_v=t(i,l);     % global dof-number of vertex i of actual element
      for j=1:3           % loop over local dofs (index of ansatz function v) 
        jglob_v=t(j,l);   % global dof-number of vertex j of actual element 

        % assemble mass matrix M
        dM=Jdet*weight(k)*sh_ref(j,k)*sh_ref(i,k);
        M(iglob_v,jglob_v)=M(iglob_v,jglob_v)+dM;

      end
    end                                   
  end
end

return











