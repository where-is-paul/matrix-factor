function [A,M]=assemble_maxwell(p,t,te,edge_num)
%

% This function assembles the global stiffness and mass matrices
% for the curl curl + I operator using conforming Nedelec_1-elements 
% of lowest order. NOTE: The boundary conditions are not taken
% into acccount.

disp('assembling stiffness and mass matrices for Nedelec elements....')

% get point and element data 
%
el_num =length(t);                         % number of elements
vertex_num = length(p);                    % number of vertices
x=p(1,:);                                  % global x-coordinates of vertices
y=p(2,:);                                  % global y-coordinates of vertices

A=sparse(edge_num,edge_num);
M=sparse(edge_num,edge_num);

% quadrature points on master-element
%
[xnod,ynod,weight]=quad5;
quad_num = length(weight);

% get Nedelec reference shape functions
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
 
  J=[x1-x0,x2-x0;y1-y0,y2-y0];  % Jacobian of affine element-map      
  trans=[x0;y0];                % translation vector of affine element-map
  JinvT=transpose(inv(J));      % transpose of inverse  
  Jdet=det(J);                  % determinant of J
  
  % element curlcurl-matrix: Ak = 2/Jdet * ones(3)

  Ak = 2/Jdet;
  
  for k=1:quad_num                % loop over quadrature points
    dx=Jdet*weight(k);            % volume element  

    [quad]=J*[xnod(k);ynod(k)]+trans;   % physical point
    xquad=quad(1);  
    yquad=quad(2); 
    
    for i=1:3  % loop over local dofs (index of test function v)
                       
      iglob_e=te(i,l);   % global dof-number of edge i of actual element 
       
      %%%% Piola-transformation for i-th Nedelec shape-function
      %%%% evaluated at quadrature point (x(k),y(k)) %%%%%
     
      shi=JinvT*sh_ref(:,i,k);      
   
      for j=1:3   % loop over local dofs (index of trial functions u and p)

        jglob_e=te(j,l);    % global dof-number of edge j of actual element 

        %%%% Piola-transformation for j-th Nedelec shape-function
        %%%% evaluated at quadrature point (x(k),y(k)) %%%%
        
        shj=JinvT*sh_ref(:,j,k); 
         
        % assemble stiffness-matrix A (pay attention to orientation of edge!)%
        %
        A(iglob_e,jglob_e)=A(iglob_e,jglob_e)+te(i+3,l)*te(j+3,l)*Ak/quad_num;
        
        % assemble mass-matrix 
        %
        M(iglob_e,jglob_e)=M(iglob_e,jglob_e)+te(i+3,l)*te(j+3,l)*transpose(shi)*shj*dx;
 
     end   
    end
  end                                   
end

return

