function [B]=assemble_coupling(p,t,te,edge_num)

% This function assembles the coupling-matrix B using
% continuous Nedelec_1-elements of lowest order 
% (linear in each component) for u and linear nodal elements for p. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get point and element data 

el_num =length(t);                         % number of elements
vertex_num = length(p);                    % number of vertices
x=p(1,:);                                  % global x-coordinates of vertices
y=p(2,:);                                  % global y-coordinates of vertices


B=sparse(edge_num,vertex_num);

% quadrature points on master-element
%
[xnod,ynod,weight]=quad5;
% number of quadrature points
%
quad_num = length(weight);

% get Nedelec reference shape functions
%
sh_ref=Ned1_0(xnod,ynod);

% gradients of nodal reference shape functions:

ref_nodal_grads = [[-1;-1],[1;0],[0;1]];


for l=1:el_num                            % loop over elements
 
  % get global coordinates of vertices 1, 2, 3 of current element
  x0=x(t(1,l));   
  x1=x(t(2,l));
  x2=x(t(3,l));
  y0=y(t(1,l));
  y1=y(t(2,l));
  y2=y(t(3,l));
 
  J=[x1-x0,x2-x0;y1-y0,y2-y0];        % Jacobian       
  trans=[x0;y0];                      % translation vector 
  JinvT=transpose(inv(J));            % transpose of inverse
  Jdet=det(J);                        % determinant of J
  
  % Transformed gradients of nodal shape functions
  nodal_grads = JinvT*ref_nodal_grads;
 
  for k=1:quad_num                % loop over quadrature points
    
    dx=Jdet*weight(k);            % volume element  

    for i=1:3                     % loop over edge shape functions 
      
        %%%% Piola-transformation for i-th Nedelec shape-function
        %%%% evaluated at quadrature point (x(k),y(k)) %%%%
       
        shi=JinvT*sh_ref(:,i,k);  
        
        iglob_e=te(i,l);                  % global number of edge
      
        for j=1:3                         % loop over test dofs
	  jglob_v=t(j,l);  
             
  B(iglob_e,jglob_v)=B(iglob_e,jglob_v) + te(i+3,l)*shi'*nodal_grads(:,j)*dx;
       
        end   
    end
  end                                   
end


return



