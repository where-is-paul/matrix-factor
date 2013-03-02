assemble_coupling.m                                                                                 0000700 0116504 0001751 00000004235 11757300770 014265  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [B]=assemble_coupling(p,t,te,edge_num)

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



                                                                                                                                                                                                                                                                                                                                                                   assemble_laplacian.m                                                                                0000700 0116504 0001751 00000003674 11757300770 014377  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [S]=assemble_scalar_lapacian(p,t)
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



                                                                    assemble_load.m                                                                                     0000700 0116504 0001751 00000004223 11757300770 013361  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [L]=assemble_load(p,t,te,edge_num)
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








                                                                                                                                                                                                                                                                                                                                                                             assemble_maxwell.m                                                                                  0000700 0116504 0001751 00000005237 11757300770 014121  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [A,M]=assemble_maxwell(p,t,te,edge_num)
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

                                                                                                                                                                                                                                                                                                                                                                 assemble_scalar_laplacian.m                                                                         0000700 0116504 0001751 00000004012 11757300770 015707  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [A]=assemble_scalar_laplacian(p,t)

%
% assembles a continuous P1 stiffness matrix: A
%                                load vector: L
% insert right hand side below


%%%%%%%%%%%%%%%%%%%%% list of points %%%%%%%%%%%%%%%%%%

Nnod=length(p)                          % number of nodes
Nel =length(t)                          % number of elements
x=p(1,:);    
y=p(2,:);


L=zeros(Nnod,1);
A=sparse(Nnod,Nnod);

%%%%%%%%%%%%%%%%%%%%% quadrature points %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xnod,ynod,weight]=quad5;

%%%%%%%%%%%%%%%%%%%%% shape functions and their derivatives %%%%%%%%%

sh=shap1(xnod,ynod);
[dxsh,dysh]=dshap1(xnod,ynod);

display('assembling Laplacian...');

for l=1:Nel                            % loop over elements
 
  % elemental data
  x0=x(t(1,l));   
  x1=x(t(2,l));
  x2=x(t(3,l));
  y0=y(t(1,l));
  y1=y(t(2,l));
  y2=y(t(3,l));
 
  B=[x1-x0,x2-x0;y1-y0,y2-y0];        % Jacobian of elemental map      
  trans=[x0;y0];                      % translation vector 
  Binv=inv(B);                        % inverse
  Jdet=det(B);                        % determinant of B

  for k=1:length(weight)                % loop over quadrature points
    
    dv=Jdet*weight(k);                  % volume element  

    [pphys]=B*[xnod(k);ynod(k)]+trans;  % physical point
    xp=pphys(1);  
    yp=pphys(2); 
 
    for j=1:3                           % loop over dofs
                       
      jj=t(j,l);                        % global dof 
                                        % gradient at integration point 
      gradj=transpose(Binv)*[dxsh(j,k);dysh(j,k)];
   
      for i=1:3                         % loop over dofs 
      
        ii=t(i,l);                      % global dof
                                        % gradient at integration point
        gradi=transpose(Binv)*[dxsh(i,k);dysh(i,k)];
      
        % assemble stiffness matrix

        A(jj,ii)=A(jj,ii)+gradj(1)*gradi(1)*dv+gradj(2)*gradi(2)*dv;

      end   
    end   
      
  end                                   
end



































                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      assemble_scalar_load.m                                                                              0000700 0116504 0001751 00000004270 11757300770 014710  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [L]=assemble_scalar_load(p,t,te,node_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assembles the load vector for the Laplace discretization with
% piecewise linear. Boundary conditions are not taken
% into account.
% NOTE: Data has to be adjusted in data_scalar.m and exact_scalar.m
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








                                                                                                                                                                                                                                                                                                                                        assemble_scalar_mass.m                                                                              0000700 0116504 0001751 00000003121 11757300770 014726  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [M]=assemble_scalar_mass(p,t)

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











                                                                                                                                                                                                                                                                                                                                                                                                                                               convergence_mixed.m                                                                                 0000700 0116504 0001751 00000006177 11757300770 014265  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program tests the convergence properties of our mixed implementation
%
% NOTE: The data has to be adjusted in data_m.m and the exact solution is
% in exact_m.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('adjust the data and exact solution in data.m and exact.m')
  

ref=4;                      % number of global refinements
clear L2_error_u   
clear curl_error_u 
clear L2_error_p   
clear H1_error_p   
clear l_num        
clear DoFs        

% initialize a mesh
[p,e,t]=initmesh('squareg','hmax',1.5);

for i=1:ref

  i

  [p,e,t]=refinemesh('squareg',p,e,t,'regular');
    
  % find and number the edges
  [te, be, edge_num]=edges(t,e);

  % identify interior/boundary edges
  [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num);

  % identify interior/boundary nodes 
  [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e);
 
  % dimensions
  %
  int_edge_num=length(list_of_int_edges);
  bnd_edge_num=length(list_of_bnd_edges);
  int_node_num =length(list_of_int_nodes);
  bnd_node_num = length(list_of_bnd_nodes);
  edge_num=int_edge_num+bnd_edge_num;
  node_num=int_node_num+bnd_node_num;

  % permutations of interior/boundary dofs

  clear edge_index_reverse; 
  edge_index=[list_of_int_edges,list_of_bnd_edges];
  edge_index_reverse(edge_index(:))=1:edge_num;

  clear node_index_reverse;
  node_index=[list_of_int_nodes,list_of_bnd_nodes];
  node_index_reverse(node_index(:))=1:node_num;

  % get our matrices
  [A11,M11,L1,S11,Q11,B11,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Mixed System
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  K = [A11 B11; B11' sparse(int_node_num,int_node_num)];
  F = [ L1; zeros(int_node_num,1)];

  %%% solution of the KKT system

  cond(full(A11+M11)),pause
  [sol] = K\F;



  %%%% electric field and potential (only interior)


  u_i = sol(1:int_edge_num,1);
  p_i = sol(int_edge_num+1:length(sol), 1);

  u_per = [u_i;u_b];
  p_per = [p_i;p_b];

  % inverse permutation of interior/boundary edges
  edge_index=[list_of_int_edges,list_of_bnd_edges];
  node_index=[list_of_int_nodes,list_of_bnd_nodes];


  uh=u_per(edge_index_reverse);
  ph=p_per(node_index_reverse);



  [L2_error_u(i),curl_error_u(i)] = uh_errors(p,t,te,uh)
  [L2_error_p(i),H1_error_p(i)] = ph_errors(p,t,ph)

  % compute number of elements
  l_num(i) = length(t);
    
  % compute number of degrees of freedom in u
  DoFs(i) = length(uh);


end



figure(1)
subplot(2,1,1)
loglog(l_num,L2_error_u,'o-',l_num,curl_error_u,'*-',l_num,l_num.^(-1/2),'r')
title('L2- and H(curl)-error of uh vs. number of elements')
legend('L2-error','H(curl)-error','Order 1')
xlabel('number of elements in mesh')
ylabel('errors')
grid on

subplot(2,1,2)
loglog(l_num,L2_error_p,'o-',l_num,H1_error_p,'*-',l_num,l_num.^(-1/2),l_num,l_num.^(-1),'-.' )
title('L2- and H1-error of ph vs. number of elements')
legend('L2-error','H1-error','Order 1', 'Order 2')
xlabel('number of elements in mesh')
ylabel('errors')
grid on







                                                                                                                                                                                                                                                                                                                                                                                                 data.m                                                                                              0000700 0116504 0001751 00000001146 11757300770 011501  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                                                                                                                                                                                                                                                                                                                                                                                                                          driver.m                                                                                            0000700 0116504 0001751 00000000210 11757300770 012052  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                main
[A11,M11,L1,S11,Q11,B11,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,...
list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);
                                                                                                                                                                                                                                                                                                                                                                                        dshap1.m                                                                                            0000700 0116504 0001751 00000000612 11757300770 011745  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [dxsh,dysh]=dshap1(x,y)
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


                                                                                                                      edges.m                                                                                             0000700 0116504 0001751 00000011355 11757300770 011662  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [Te, Be, edge_ind] = edges(T,E)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Anna S., 4.10.02
%
% this function numbers the edges of a triangular mesh 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The input parameter is the element-list T, as generated by the PDE-toolbox-command 'initmesh'
% (for description type 'help initmesh' in the MATLAB command-window). In fact, we only need the first 3 rows
% of T.
%
% The output parameters are the total number of edges edge_ind, the matrix
% Te containing the edge data and the vector Be containing the indices of
% the boundary-edges.
% The structure of Te is similar to the structure of T: the 3 first rows of Te contain the global edge-numbers of
% the elements (each column of Te represents an element). Locally on an element, the edges are numbered couterclockwise, 
% starting with the one between the local vertices 1 and 2. The 3 last rows of Te endow the edges with an orientation,
% which is set to plus or minus one. 
% 
% The numbering is straight-forward and not sophisticated at all; we cannot expect the non-zero entries of the 
% global stiffness-matrix A to be disributed in a regular pattern. This can affect the condition number of A and the
% memory requirement in a negative way.
% 
%
% Accelerated version: October 02


disp('computing edge data...')
pause(0.001)        % this is just a hack to force MATLAB to display the previous command immediately



el_num = length(T);        % number of elements
edge_ind = 0;               % edge-index
Te = zeros(6,el_num);      % number of edges are 3 per element (+3 rows for edge-orientation)
Be = [];                % initialize boundary edge vector
be_num = length(E);        % number of boundary edges


% loop over all elements
for l=1:el_num
% loop over edges of element
for e=1:3
    
    % if edge is not yet numbered...
    if Te(e,l)== 0 
        % number edge
        edge_ind = edge_ind + 1;
        Te(e,l) = edge_ind;
        % give orientation +1
        Te(e+3,l) = 1;
        
        % get global number of start- and end-point
        Ps = T(e,l);
        Pe = T(mod(e,3)+1,l);
        
        %% forward-search for other element containing this edge %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        stop = 0;                % control variable to stop the numbering if other edge is found
        for ol = l+1 : el_num
            for Ps_local = 1 : 3
                if ( ~(T(Ps_local,ol)-Ps) | ~(T(Ps_local,ol)-Pe) )    % start- or end-point are in element ol
                    
                    for Pe_local = (Ps_local+1) : 3
                       if ( ~(T(Pe_local,ol)-Ps) | ~(T(Pe_local,ol)-Pe) ) 
                           % other edge between local pts Ps_local and Pe_local
                           oe_local = Ps_local + Pe_local;
                           switch oe_local
                           case 3 
                            % edge between local pts 1 and 2
                            Te(1,ol) = edge_ind;
                            Te(4,ol) = -1;
                           case 5
                            % edge between local pts 2 and 3
                            Te(2,ol) = edge_ind;
                            Te(5,ol) = -1;
                           case 4
                            % edge between local pts 3 and 1
                            Te(3,ol) = edge_ind;
                            Te(6,ol) = -1;
                           end
                           stop = 1;        % other edge is found and we can stop this nested loop
                       end  % if
                   end  % for Pe_local
               end  % if
           end % for Ps_local
           
           if(stop)         % leave outer loop, since there is at most 1 other element containing e as edge
                break    
           end
               
           end  % for ol
           
           %%% In order to have an easy access to boundary edges, we flag them here %%%
           
           for be = 1 : be_num      % loop over all boundary edges
               
               % check if current edge is a boundary edge
               if( (~(E(1,be)-Ps) | ~(E(1,be)-Pe) ) & (~(E(2,be)-Ps) | ~(E(2,be)-Pe) ) )
                   
                   %check, if this edge is on \emph{outer} boundary
                   
                    if E(6,be) | E(7,be)
                        
                        % label this edge as Dirichlet boundary edge
                        Be = [Be edge_ind];
               
                    end % if
                    break       % edge flagged as boundary edge, exit loop over boundary edges.
                end % if
                
            end % for be

        
    end % if
end % for e
end %for l


return
        
        
            







                                                                                                                                                                                                                                                                                   exact.m                                                                                             0000700 0116504 0001751 00000001510 11757300770 011667  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [u, curlu, p, gradp] = exact(x,y)
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






                                                                                                                                                                                        genMat_backup.m                                                                                     0000700 0116504 0001751 00000003275 11757300770 013335  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:4,
    i
    if i==1 
        % initialize a mesh
        [p,e,t]=initmesh('squareg','hmax',1.5);
    else
        [p,e,t]=refinemesh('squareg',p,e,t,'regular');
    end
    %p=jigglemesh(p,e,t);
    %pdeplot(p,e,t);
    %pause
    
    % find and number the edges
    [te, be, edge_num]=edges(t,e);

    % identify interior/boundary edges
    [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num);

    % identify interior/boundary nodes 
    [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e);

    % dimensions
    %
    int_edge_num=length(list_of_int_edges);
    bnd_edge_num=length(list_of_bnd_edges);
    int_node_num =length(list_of_int_nodes);
    bnd_node_num = length(list_of_bnd_nodes);
    edge_num=int_edge_num+bnd_edge_num;
    node_num=int_node_num+bnd_node_num;

    % permutations of interior/boundary dofs

    clear edge_index_reverse; 
    edge_index=[list_of_int_edges,list_of_bnd_edges];
    edge_index_reverse(edge_index(:))=1:edge_num;

    clear node_index_reverse;
    node_index=[list_of_int_nodes,list_of_bnd_nodes];
    node_index_reverse(node_index(:))=1:node_num;

    % get our matrices
    %

    [A11,M11,L1,S11,Q11,B11,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);


    B11=B11';

    [m,n]=size(B11);
    sys_size=m+n, nel=size(t,2)
    
    cond(full(A11+M11))
    
    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    %AW=A+B'*(W\B);
    %AI=A+B'*B;
    
    %save(['lshape_matrices',num2str(i),'.mat']);   
end








                                                                                                                                                                                                                                                                                                                                   genMatLshape.m                                                                                      0000700 0116504 0001751 00000002547 11757300770 013146  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                
disp('adjust the data and exact solution in data.m and exact.m')
  

% initialize a mesh
%[p,e,t]=initmesh('lshapeg','hmax',.1);

%pdeplot(p,e,t),pause

for ii=1:5

  ii

  %[p,e,t]=adaptmesh('lshapeg',p,e,t,'regular');
  load(['Lshape_mesh',num2str(ii),'.mat']); 
  
  % find and number the edges
  [te, be, edge_num]=edges(t,e);

  % identify interior/boundary edges
  [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num);

  % identify interior/boundary nodes 
  [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e);
 
  % dimensions
  %
  int_edge_num=length(list_of_int_edges);
  bnd_edge_num=length(list_of_bnd_edges);
  int_node_num =length(list_of_int_nodes);
  bnd_node_num = length(list_of_bnd_nodes);
  edge_num=int_edge_num+bnd_edge_num;
  node_num=int_node_num+bnd_node_num;

  % permutations of interior/boundary dofs

  clear edge_index_reverse; 
  edge_index=[list_of_int_edges,list_of_bnd_edges];
  edge_index_reverse(edge_index(:))=1:edge_num;

  clear node_index_reverse;
  node_index=[list_of_int_nodes,list_of_bnd_nodes];
  node_index_reverse(node_index(:))=1:node_num;

  % get our matrices
  [A,M,f,L,W,B,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);
  B=B';
  [m,n]=size(B);
  sys_size=m+n, nel=size(t,2)
      
  save(['Lshape_matrices',num2str(ii),'.mat']);   

end   








                                                                                                                                                         genMat.m                                                                                            0000700 0116504 0001751 00000002443 11757300770 012004  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                
disp('adjust the data and exact solution in data.m and exact.m')
  

% initialize a mesh
[p,e,t]=initmesh('squareg','hmax',1);

for i=1:2

  i

  [p,e,t]=refinemesh('squareg',p,e,t,'regular');
  
  % find and number the edges
  [te, be, edge_num]=edges(t,e);

  % identify interior/boundary edges
  [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num);

  % identify interior/boundary nodes 
  [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e);
 
  % dimensions
  %
  int_edge_num=length(list_of_int_edges);
  bnd_edge_num=length(list_of_bnd_edges);
  int_node_num =length(list_of_int_nodes);
  bnd_node_num = length(list_of_bnd_nodes);
  edge_num=int_edge_num+bnd_edge_num;
  node_num=int_node_num+bnd_node_num;

  % permutations of interior/boundary dofs

  clear edge_index_reverse; 
  edge_index=[list_of_int_edges,list_of_bnd_edges];
  edge_index_reverse(edge_index(:))=1:edge_num;

  clear node_index_reverse;
  node_index=[list_of_int_nodes,list_of_bnd_nodes];
  node_index_reverse(node_index(:))=1:node_num;

  % get our matrices
  [A,M,f,L,W,B,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);
  B=B';
  [m,n]=size(B);
  sys_size=m+n, nel=size(t,2)
      
  save(['matrices_homogeneous',num2str(i),'.mat']);   

end   








                                                                                                                                                                                                                             genMatSquare.m                                                                                      0000700 0116504 0001751 00000002427 11757300770 013167  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                
disp('adjust the data and exact solution in data.m and exact.m')
  

% initialize a mesh
[p,e,t]=initmesh('squareg','hmax',1);

for i=1:8

  i

  [p,e,t]=refinemesh('squareg',p,e,t,'regular');
  
  % find and number the edges
  [te, be, edge_num]=edges(t,e);

  % identify interior/boundary edges
  [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num);

  % identify interior/boundary nodes 
  [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e);
 
  % dimensions
  %
  int_edge_num=length(list_of_int_edges);
  bnd_edge_num=length(list_of_bnd_edges);
  int_node_num =length(list_of_int_nodes);
  bnd_node_num = length(list_of_bnd_nodes);
  edge_num=int_edge_num+bnd_edge_num;
  node_num=int_node_num+bnd_node_num;

  % permutations of interior/boundary dofs

  clear edge_index_reverse; 
  edge_index=[list_of_int_edges,list_of_bnd_edges];
  edge_index_reverse(edge_index(:))=1:edge_num;

  clear node_index_reverse;
  node_index=[list_of_int_nodes,list_of_bnd_nodes];
  node_index_reverse(node_index(:))=1:node_num;

  % get our matrices
  [A,M,f,L,W,B,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);
  B=B';
  [m,n]=size(B);
  sys_size=m+n, nel=size(t,2)
      
  save(['matrices',num2str(i),'.mat']);   

end   








                                                                                                                                                                                                                                         genPC2_Lshape.m                                                                                     0000700 0116504 0001751 00000001262 11757300770 013141  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
for i=1:5,
    i
    %disp('Loading matrices...');
    load(['Lshape_matrices',num2str(i),'.mat']); 
    [m,n]=size(B);
    disp('[m; n; n+m]');
    [m; n; m+n]
    Q=[A-k^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];
    %disp('Setting up preconditioners...');
    P =[A+M  sparse(n,m); sparse(m,n) L];

    %disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); 
    %flag, iter,relres,resvec
    iter
    %save(['Lshape_preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                              genPC2.m                                                                                            0000700 0116504 0001751 00000001233 11757300770 011643  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=1/8;
for i=3:3,
    i
    disp('Loading matrices...');
    load(['matrices',num2str(i),'.mat']); 
    [m,n]=size(B);
    disp('[m; n; n+m]');
    [m; n; m+n]
    Q=[A-k^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];
    disp('Setting up preconditioners...');
    P =[A+M  sparse(n,m); sparse(m,n) L];

    disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); 
    flag, iter,relres,resvec
 
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                                                     genPC3.m                                                                                            0000700 0116504 0001751 00000001270 11757300770 011645  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=0;
for i=2:2,
    i
    disp('Loading matrices...');
    load(['matrices',num2str(i),'.mat']); 
    [m,n]=size(B);
    disp('[m; n; n+m]');
    [m; n; m+n]
    Q=[A-k^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];
    disp('Setting up preconditioners...');
    P =[A+M  sparse(n,m); sparse(m,n) L];

    d=eig(full(Q),full(P));
    %disp('Solving with MINRES...')
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); 
    %flag, iter,relres,resvec
 
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                        genPC_homogeneous_inexact.m                                                                         0000700 0116504 0001751 00000004226 11757300770 015711  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=0;

for i=1:3
    i
    disp('Loading matrices...');
    load(['matrices_homogeneous',num2str(i),'.mat']); 
    
    [m,n]=size(B);
    Q=[A-kk^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];

    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    
    %disp('Running PCG');
    %[xcg,flagcg,relrescg,itercg,resveccg] = mypcg_inexact(A-kk^2*M,f,1e-10,1000,A+(1-kk^2)*M,[],[],L,B); flagcg,itercg
    
    
    %AW=A+B'*(W\B);
    
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AL,f,1e-12,100,A+M); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AW,f,1e-12,100,A+nel*B'*B); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AI,f,1e-10,100,A+); ITER,FLAG
    %pause
    %disp('Setting up preconditioners...');
    P =[A+(1-kk^2)*M  sparse(n,m); sparse(m,n) L];
    %P2 =[A+M  sparse(n,m); sparse(m,n) W];
    %PL=[AL-kk^2*M sparse(n,m); sparse(m,n) L];
    %PW=[AW-kk^2*M sparse(n,m); sparse(m,n) W];
    %I=speye(m);
    %delta=1/nel;
    %PI=[A-kk^2*M+1/delta*B'*B sparse(n,m);sparse(m,n) delta*I];
    
    %disp('Computing eigenvalues...');
    %eig_P  = sort(real(eig(full(Q),full(P)))),pause
   %eig_PL = sort(real(eig(full(Q),full(PL)))),pause
    %eig_PW = sort(real(eig(full(Q),full(PW)))),pause
    %eig_P2 = sort(real(eig(full(Q),full(P2)))),pause
    %eig_PI = sort(real(eig(full(Q),full(PI))));
    %disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); flag, iter
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,PI); flag, iter
    %[x2,flag2,relres2,iter2,resvec2,resveccg2] = minres(Q,rhs,1e-11,500,P2); flag2, iter2
    
    %[xL,flagL,relresL,iterL,resvecL,resveccgL] = minres(Q,rhs,1e-10,500,PL); flagL, iterL
    %[xW,flagW,relresW,iterW,resvecW,resveccgW] = minres(Q,rhs,1e-10,500,PW); flagW, iterW
    
    %pause
    %DUMMY = sparse(m+n,m+n);
    %tic,[xW2,flagW2,relresW2,iterW2,resvecW2,resveccgW2] = myminres(Q,rhs,1e-10,500,DUMMY,[],[],L,A,B); flagW2,iterW2,toc 
    
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                                                          genPC_homogeneous.m                                                                                 0000700 0116504 0001751 00000004213 11757300770 014172  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=0

for i=1:7
    i
    disp('Loading matrices...');
    load(['matrices_homogeneous',num2str(i),'.mat']); 
    
    [m,n]=size(B);
    Q=[A-kk^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];

    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    
    disp('Running PCG');
    [xcg,flagcg,relrescg,itercg,resveccg] = mypcg(A-kk^2*M,f,1e-10,1000,A+(1-kk^2)*M,[],[],L,B); flagcg,itercg
    
    
    %AW=A+B'*(W\B);
    
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AL,f,1e-12,100,A+M); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AW,f,1e-12,100,A+nel*B'*B); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AI,f,1e-10,100,A+); ITER,FLAG
    %pause
    %disp('Setting up preconditioners...');
    P =[A+(1-kk^2)*M  sparse(n,m); sparse(m,n) L];
    %P2 =[A+M  sparse(n,m); sparse(m,n) W];
    %PL=[AL-kk^2*M sparse(n,m); sparse(m,n) L];
    %PW=[AW-kk^2*M sparse(n,m); sparse(m,n) W];
    %I=speye(m);
    %delta=1/nel;
    %PI=[A-kk^2*M+1/delta*B'*B sparse(n,m);sparse(m,n) delta*I];
    
    %disp('Computing eigenvalues...');
    %eig_P  = sort(real(eig(full(Q),full(P)))),pause
   %eig_PL = sort(real(eig(full(Q),full(PL)))),pause
    %eig_PW = sort(real(eig(full(Q),full(PW)))),pause
    %eig_P2 = sort(real(eig(full(Q),full(P2)))),pause
    %eig_PI = sort(real(eig(full(Q),full(PI))));
    %disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); flag, iter
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,PI); flag, iter
    %[x2,flag2,relres2,iter2,resvec2,resveccg2] = minres(Q,rhs,1e-11,500,P2); flag2, iter2
    
    %[xL,flagL,relresL,iterL,resvecL,resveccgL] = minres(Q,rhs,1e-10,500,PL); flagL, iterL
    %[xW,flagW,relresW,iterW,resvecW,resveccgW] = minres(Q,rhs,1e-10,500,PW); flagW, iterW
    
    %pause
    %DUMMY = sparse(m+n,m+n);
    %tic,[xW2,flagW2,relresW2,iterW2,resvecW2,resveccgW2] = myminres(Q,rhs,1e-10,500,DUMMY,[],[],L,A,B); flagW2,iterW2,toc 
    
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                                                                     genPC_Lshape.m                                                                                      0000700 0116504 0001751 00000004060 11757300770 013056  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=1/8;

for ii=1:5
    ii
    disp('Loading matrices...');
    load(['Lshape_matrices',num2str(ii),'.mat']); 
    [m,n]=size(B);
    Q=[A-kk^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];

    
    disp('Running PCG');
    [xcg,flagcg,relrescg,itercg,resveccg] = mypcg(A-kk^2*M,f,1e-10,1000,A+(1-kk^2)*M,[],[],L,B); flagcg,itercg
    
    
    
    %condest(L),pause
    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    
    %condest(AL)
    %condest(AW),pause
    
    %AW=A+B'*(W\B);
    %W=normest(A)/normest(B'*B)*speye(size(B,1));
    %nAB=norm(A,1)/(norm(B,1)^2);
    %AW=A+B'*(W\B);
    
    
    
    %disp('Setting up preconditioners...');
    %P =[A+(1-kk^2)*M  sparse(n,m); sparse(m,n) L];
    %PL=[AL-k^2*M sparse(n,m); sparse(m,n) L];
    %PW=[AW-k^2*M sparse(n,m); sparse(m,n) W];
    
    %I=speye(m);
    %delta=1/nel;
    %PI=[A-kk^2*M+B'*B sparse(n,m);sparse(m,n) I];
    
    %[xcg,flagcg,relrescg,itercg,resveccg] = PCG(AL,f,1e-10,1000,A+(1-kk^2)*M); flagcg,itercg
    
    %condest(PL),pause
    
    %sort(eig(full(Q),full(PW))),pause
    
    %disp('Computing eigenvalues...');
    %eig_P  = sort(real(eig(full(Q),full(P))));
    %eig_PL = sort(real(eig(full(Q),full(PL))));
    %eig_PW = sort(real(eig(full(Q),full(PW))));
    %eig_PI = sort(real(eig(full(Q),full(PI))));
    %disp('Solving with MINRES...')
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); flag, iter
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,PI); flag, iter
    %[xL,flagL,relresL,iterL,resvecL,resveccgL] = minres(Q,rhs,1e-6,500,PL); flagL, iterL
    %tic,[xW,flagW,relresW,iterW,resvecW,resveccgW] = minres(Q,rhs,1e-10,500,PW); flagW, iterW,toc
    
    %DUMMY = sparse(m+n,m+n);
    %tic,[xL2,flagL2,relresL2,iterL2,resvecL2,resveccgL2] = myminres(Q,rhs,1e-8,500,DUMMY,[],[],L,A,B); flagL2,iterL2,toc 
    
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                                                                                                                                                                genPC.m                                                                                             0000700 0116504 0001751 00000004202 11757300770 011560  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kk=1/8;

for i=1:6
    i
    disp('Loading matrices...');
    load(['matrices',num2str(i),'.mat']); 
    
    [m,n]=size(B);
    Q=[A-kk^2*M B'; B sparse(m,m)];
    rhs=[f; zeros(m,1)];

    %disp('Computing the augmented matrices...')
    %AL=A+B'*(L\B);
    
    disp('Running PCG');
    [xcg,flagcg,relrescg,itercg,resveccg] = mypcg(A-kk^2*M,f,1e-10,1000,A+(1-kk^2)*M,[],[],L,B); flagcg,itercg
    
    
    %AW=A+B'*(W\B);
    
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AL,f,1e-12,100,A+M); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AW,f,1e-12,100,A+nel*B'*B); ITER,FLAG
    %[X,FLAG,RELRES,ITER,RESVEC] = pcg(AI,f,1e-10,100,A+); ITER,FLAG
    %pause
    %disp('Setting up preconditioners...');
    P =[A+(1-kk^2)*M  sparse(n,m); sparse(m,n) L];
    %P2 =[A+M  sparse(n,m); sparse(m,n) W];
    %PL=[AL-kk^2*M sparse(n,m); sparse(m,n) L];
    %PW=[AW-kk^2*M sparse(n,m); sparse(m,n) W];
    %I=speye(m);
    %delta=1/nel;
    %PI=[A-kk^2*M+1/delta*B'*B sparse(n,m);sparse(m,n) delta*I];
    
    %disp('Computing eigenvalues...');
    %eig_P  = sort(real(eig(full(Q),full(P)))),pause
   %eig_PL = sort(real(eig(full(Q),full(PL)))),pause
    %eig_PW = sort(real(eig(full(Q),full(PW)))),pause
    %eig_P2 = sort(real(eig(full(Q),full(P2)))),pause
    %eig_PI = sort(real(eig(full(Q),full(PI))));
    %disp('Solving with MINRES...')
    [x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,P); flag, iter
    %[x,flag,relres,iter,resvec,resveccg] = minres(Q,rhs,1e-10,500,PI); flag, iter
    %[x2,flag2,relres2,iter2,resvec2,resveccg2] = minres(Q,rhs,1e-11,500,P2); flag2, iter2
    
    %[xL,flagL,relresL,iterL,resvecL,resveccgL] = minres(Q,rhs,1e-10,500,PL); flagL, iterL
    %[xW,flagW,relresW,iterW,resvecW,resveccgW] = minres(Q,rhs,1e-10,500,PW); flagW, iterW
    
    %pause
    %DUMMY = sparse(m+n,m+n);
    %tic,[xW2,flagW2,relresW2,iterW2,resvecW2,resveccgW2] = myminres(Q,rhs,1e-10,500,DUMMY,[],[],L,A,B); flagW2,iterW2,toc 
    
    %save(['preconditioners',num2str(i),'.mat']); 
    
end








                                                                                                                                                                                                                                                                                                                                                                                              iterchk.m                                                                                           0000700 0116504 0001751 00000002472 11757300770 012224  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                
function [atype,afun,afcnstr] = iterchk(A)
%ITERCHK  Checks arguments to iterative methods.
%   [ATYPE,AFUN,AFCNSTR] = ITERCHK(A) returns the following:
%   ATYPE is either 'matrix', 'function', 'expression' or 'inline object'.
%   AFUN is the function name or inline object.
%   AFUN is '' if ATYPE is 'matrix'.
%   AFCNSTR is the function name if ATYPE is 'function'.
%   AFCNSTR is the formula of the function if ATYPE is 'expression' or
%   'inline object'.  AFCNSTR is '' if ATYPE is 'matrix'.
%
%   See also BICG, BICGSTAB, CGS, GMRES, PCG, QMR.

%   Penny Anderson, 1998.
%   Copyright 1984-2003 The MathWorks, Inc. 
%   $Revision: 1.8.4.1 $ $Date: 2003/05/19 11:16:35 $


[afun,afunmsg] = fcnchk(A);
if isempty(afunmsg)
   if isa(afun,'inline')      
      if isa(A,'inline')
         atype = 'inline object';
      else
         atype = 'expression';
      end
      afcnstr = formula(afun);
   else % both function_handles @fun and function names 'fun'
      atype = 'function';
      if isa(A,'function_handle')
          afcnstr = func2str(A);
      else
          afcnstr = A;
      end
   end
elseif isa(A,'double')
   atype = 'matrix';
   afcnstr = '';
else
   error('MATLAB:iterchk:InvalidInput',...
         'Argument must be a matrix, a function, or an inline object.');
end

                                                                                                                                                                                                      loadMat.m                                                                                           0000700 0116504 0001751 00000001351 11757300770 012147  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
k=0;
gamma=1;
for i=1:2,
    load(['matrices',num2str(i),'.mat']); 
    [m,n]=size(B);
    F=A-k^2*M;
    Q=[F B'; B sparse(m,m)];
    P =[A+gamma*M  sparse(n,m); sparse(m,n) L];
    PL=[F+B'*(L\B) sparse(n,m); sparse(m,n) L];
    PW=[F+B'*(W\B) sparse(n,m); sparse(m,n) W];
    PI=[F+B'*B sparse(n,m); sparse(m,n) speye(m)];
    eig_P  = sort(real(eig(full(Q),full(P))));
    eig_PL = sort(real(eig(full(Q),full(PL))));
    eig_PW = sort(real(eig(full(Q),full(PW))));
    eig_PI = sort(real(eig(full(Q),full(PI))));
    save(['preconditioners',num2str(i),'.mat']); 
end








                                                                                                                                                                                                                                                                                       main.m                                                                                              0000700 0116504 0001751 00000003552 11757300770 011517  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% initialize a mesh
[p,e,t]=initmesh('lshapeg','hmax',0.1);

% find and number the edges
[te, be, edge_num]=edges(t,e);

% identify interior/boundary edges
[list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num);

% identify interior/boundary nodes 
[list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e);

% dimensions
%
int_edge_num=length(list_of_int_edges);
bnd_edge_num=length(list_of_bnd_edges);
int_node_num =length(list_of_int_nodes);
bnd_node_num = length(list_of_bnd_nodes);
edge_num=int_edge_num+bnd_edge_num;
node_num=int_node_num+bnd_node_num;

% permutations of interior/boundary dofs

clear edge_index_reverse; 
edge_index=[list_of_int_edges,list_of_bnd_edges];
edge_index_reverse(edge_index(:))=1:edge_num;

clear node_index_reverse;
node_index=[list_of_int_nodes,list_of_bnd_nodes];
node_index_reverse(node_index(:))=1:node_num;

% get our matrices
[A11,M11,L1,S11,Q11,B11,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Mixed System
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = [A11 B11; B11' sparse(int_node_num,int_node_num)];
F = [ L1; zeros(int_node_num,1)];

%%% solution of the KKT system

[sol] = K\F;



%%%% electric field and potential (only interior)


u_i = sol(1:int_edge_num,1);
p_i = sol(int_edge_num+1:length(sol), 1);

u_per = [u_i;u_b];
p_per = [p_i;p_b];

% inverse permutation of interior/boundary edges
edge_index=[list_of_int_edges,list_of_bnd_edges];
node_index=[list_of_int_nodes,list_of_bnd_nodes];


uh=u_per(edge_index_reverse);
ph=p_per(node_index_reverse);


[l2_error_u,curl_error_u]=uh_errors(p,t,te,uh)
[l2_error_p,h1_error_p]=ph_errors(p,t,ph)

return









                                                                                                                                                      myadaptmesh.m                                                                                       0000700 0116504 0001751 00000024351 11757300770 013107  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                

function [u,p,e,t]=adaptmesh(g,b,c,a,f,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24)
%ADAPTMESH Adaptive mesh generation and PDE solution.
%
%       [U,P,E,T]=ADAPTMESH(G,B,C,A,F,P1,V1,...) performs adaptive mesh
%       generation and PDE solution.  The large number of possible input
%       option is handled using property-value pair arguments.  The first
%       five arguments G, B, C, A, and F are not optional.
%
%       The function produces a solution U to the elliptic scalar PDE problem
%       -div(c*grad(u))+a*u=f where the problem geometry and boundary
%       conditions given by G and B.
%
%       The solution u is represented as the MATLAB column vector U.
%       See ASSEMPDE for details.
%
%       G describes the geometry of the PDE problem. G can
%       either be a Decomposed Geometry Matrix or the name of Geometry
%       M-file. See either DECSG or PDEGEOM for details.
%
%       B describes the boundary conditions of the PDE problem.  B
%       can either be a Boundary Condition Matrix or the name of Boundary
%       M-file. See PDEBOUND for details.
%
%       The adapted triangular mesh of the PDE problem is given by the triangle
%       data P, E, and T.  See either INITMESH or PDEGEOM for details.
%
%       The coefficients C, A and F of the PDE problem can
%       be given in a wide variety of ways.  See ASSEMPDE for details.
%
%       Valid property/value pairs include
%
%       Prop.   Value/{Default}         Description
%       ----------------------------------------------------------------------
%       Maxt    Positive integer {Inf}  Maximum number of new triangles
%       Ngen    Positive integer {10}   Maximum number of triangle generations
%       Mesh    P1, E1, T1              Initial mesh
%       Tripick {pdeadworst}|pdeadgsc   Triangle selection method
%       Par     Numeric {0.5}           Function parameter
%       Rmethod {longest}|regular       Triangle refinement method
%       Nonlin  on | off                Use nonlinear solver
%       Toln    numeric {1e-3}          Nonlinear tolerance
%       Init    string|numeric          Nonlinear initial solution value
%       Jac     {fixed}|lumped|full     Nonlinear solver Jacobian calculation
%       Norm    Numeric {Inf}           Nonlinear solver residual norm
%
%       Par is passed to the tripick function. Normally it is used as
%       tolerance of how well the solution fits the equation. No more than
%       Ngen successive refinements are attempted. Refinement is also
%       stopped when the number of triangles in the mesh exceeds the
%       Maxt.
%
%       P1, E1, and T1 are the input mesh data. This triangle mesh is used
%       as a starting mesh for the adaptive algorithm. If no initial mesh
%       is provided, the result of a call to INITMESH with no options is
%       used as initial mesh.
%
%       The triangle pick method is a user-definable triangle selection
%       method.  Given the error estimate computed by the function PDEJMPS,
%       the triangle pick method selects the triangles to be refined
%       in the next triangle generation. The function is called using the
%       arguments P, T, CC, AA, FF, U, ERRF, and PAR.  P and T represent the
%       current generation of triangles, CC, AA, FF are the current
%       coefficients for the PDE problem, expanded to triangle midpoints,
%       U is the current solution, ERRF is the computed error estimate,
%       and PAR, the function parameter, given to ADAPTMESH as optional
%       argument. The matrices CC, AA, FF, and ERRF all have NT columns,
%       where NT is the current number of triangles.
%       The number of rows in CC, AA, and FF are exactly the same as the
%       input arguments C, A, and F. ERRF has one row for each equation
%       in the system. There are two standard triangle selection methods
%       in the PDE Toolbox - PDEADWORST and PDEADGSC.
%       PDEADWORST selects triangles where ERRF exceeds a fraction
%       (default: 0.5) of the worst value. PDEADGSC selects triangles
%       using a relative tolerance criterion.
%
%       The refinement method is either 'longest' or 'regular'.
%       See REFINEMESH for details.
%
%       Also nonlinear PDE problems can be solved by the adaptive algorithms.
%       For nonlinear PDE problems, the 'Nonlin' parameter must be set
%       to 'on'. The nonlinear tolerance Toln and nonlinear initial value
%       U0 are passed to the nonlinear solver. See PDENONLIN for details.
%
%       See also ASSEMPDE, PDEBOUND, PDEGEOM, INITMESH, REFINEMESH, PDENONLIN

%       A. Nordmark 10-18-94, AN 01-23-95, MR 05-24-95.
%       Copyright 1994-2003 The MathWorks, Inc.
%       $Revision: 1.9.4.1 $  $Date: 2003/10/21 12:25:57 $

alfa=0.15;
beta=0.15;
mexp=1;
mesh=0;
nonl=0;
gotu=0;

% Default values
Tripick='pdeadworst';
Rmethod='longest';
Toln=1e-4;
Ngen=32;
Maxt=Inf;
Par=0.5;
Jac='fixed';
norm=Inf;

k=1;
noptarg=nargin-5;
while k<=noptarg
  Name=eval(['p' int2str(k)]);
  if ~ischar(Name)
    error('PDE:adaptmesh:ParamNotString', 'Parameter must be a string.')
  elseif size(Name,1)~=1,
    error('PDE:adaptmesh:ParamNumRowsOrEmpty', 'Parameter must be a non-empty single row string.')
  end
  Name=lower(Name);
  if strcmp(Name,'mesh')
    if noptarg-k<3
      error('PDE:adaptmesh:MeshNumValues', 'Option ''mesh'' must have three values.')
    end
    k=k+1;
    p=eval(['p' int2str(k)]);
    k=k+1;
    e=eval(['p' int2str(k)]);
    k=k+1;
    t=eval(['p' int2str(k)]);
    if ischar(p) || ischar(e) || ischar(t)
      error('PDE:adaptmesh:MeshNotNumeric', 'Mesh data must be numeric.');
    end
    mesh=1;
  elseif strcmp(Name,'tripick')
    if noptarg-k<1
      error('PDE:adaptmesh:TripickNoValue', 'Option ''tripick'' must have a value.')
    end
    k=k+1;
    Tripick=eval(['p' int2str(k)]);
    if ~ischar(Tripick)
      error('PDE:adaptmesh:TripickNotString', 'Tripick value must be a string.')
    end
  elseif strcmp(Name,'rmethod')
    if noptarg-k<1
      error('PDE:adaptmesh:RmethodNoValue', 'Option ''rmethod'' must have a value.')
    end
    k=k+1;
    Rmethod=eval(['p' int2str(k)]);
    if ~ischar(Rmethod)
      error('PDE:adaptmesh:RmethodNotString', 'Rmethod value must be a string.')
    end
  elseif strcmp(Name,'toln')
    if noptarg-k<1
      error('PDE:adaptmesh:TolnNoValue', 'Option ''toln'' must have a value.')
    end
    k=k+1;
    Toln=eval(['p' int2str(k)]);
    if ischar(Toln)
      error('PDE:adaptmesh:TolnNotNumeric', 'Toln value must be numeric.')
    end
  elseif strcmp(Name,'ngen')
    if noptarg-k<1
      error('PDE:adaptmesh:NgenNoValue', 'Option ''ngen'' must have a value.')
    end
    k=k+1;
    Ngen=eval(['p' int2str(k)]);
    if ischar(Ngen)
      error('PDE:adaptmesh:NgenNotNumeric', 'Ngen value must be numeric.')
    end
  elseif strcmp(Name,'maxt')
    if noptarg-k<1
      error('PDE:adaptmesh:MaxtNoValue', 'Option ''maxt'' must have a value.')
    end
    k=k+1;
    Maxt=eval(['p' int2str(k)]);
    if ischar(Maxt)
      error('PDE:adaptmesh:MaxtNotNumeric', 'Maxt value must be numeric.')
    end
  elseif strcmp(Name,'par')
    if noptarg-k<1
      error('PDE:adaptmesh:ParNoValue', 'Option ''par'' must have a value.')
    end
    k=k+1;
    Par=eval(['p' int2str(k)]);
    if ischar(Par)
      error('PDE:adaptmesh:ParNotNumeric', 'Par value must be numeric.')
    end
  elseif strcmp(Name,'nonlin')
    if noptarg-k<1
      error('PDE:adaptmesh:NonlinNoValue', 'Option ''nonlin'' must have a value.')
    end
    k=k+1;
    Nonlin=eval(['p' int2str(k)]);
    if ~ischar(Nonlin)
      error('PDE:adaptmesh:NonlinNotString', 'Nonlin value must be a string.')
    end
    Nonlin=lower(Nonlin);
    if strcmp(Nonlin,'on')
      nonl=1;
    elseif ~strcmp(Nonlin,'off')
      error('PDE:adaptmesh:NonlinInvalidString', 'Nonlin value must be off | {on} .')
    end
  elseif strcmp(Name,'init')
    if noptarg-k<1
      error('PDE:adaptmesh:InitNoValue', 'Option ''init'' must have a value.')
    end
    k=k+1;
    u=eval(['p' int2str(k)]);
    gotu=1;
  elseif strcmp(Name,'jac')
    if noptarg-k<1
      error('PDE:adaptmesh:JacNoValue', 'Option ''jac'' must have a value.')
    end
    k=k+1;
    Jac=eval(['p' int2str(k)]);
    if ~ischar(Jac)
      error('PDE:adaptmesh:JacNotString', 'Jac value must be a string.')
    end
    Jac=lower(deblank(Jac));
    if ~(strcmp(Jac,'fixed') || strcmp(Jac,'lumped') || strcmp(Jac,'full'))
      error('PDE:adaptmesh:JacInvalidString', 'Jac value must be {fixed} | lumped | full.')
    end
  elseif strcmp(Name,'norm')
    if noptarg-k<1
      error('PDE:adaptmesh:NormNoValue', 'Option ''norm'' must have a value.')
    end
    k=k+1;
    norm=eval(['p' int2str(k)]);
    if ischar(norm)
      error('PDE:adaptmesh:NormNotNumeric', 'Norm value must be numeric.')
    end
  else
    error('PDE:adaptmesh:InvalidOption',  ['Unknown option: ' Name])
  end
  k=k+1;
end

if ~mesh
  [p,e,t]=initmesh(g);
end
np=size(p,2);

gen=0;

jj=0;
while 1,

   
  fprintf('Number of triangles: %g\n',size(t,2))
  
  if mod(gen,8)==0
      disp('saving mesh...')
      jj=jj+1;
      save(['Lshape_mesh',num2str(jj),'.mat']);
      pdeplot(p,e,t),pause
  end
  
  
  if nonl
    if gotu
      u=pdenonlin(b,p,e,t,c,a,f,'jacobian',Jac,'U0',u,'tol',Toln,'norm',norm);
    else
      u=pdenonlin(b,p,e,t,c,a,f,'jacobian',Jac,'tol',Toln,'norm',norm);
    end
    gotu=1;
  else
    u=assempde(b,p,e,t,c,a,f);
  end

  if any(isnan(u))
    error('PDE:adaptmesh:NaNinSolution', 'Solution contains NaNs.')
  end

  % Expand values
  [cc,aa,ff]=pdetxpd(p,t,u,c,a,f);

  errf=pdejmps(p,t,cc,aa,ff,u,alfa,beta,mexp);

  i=feval(Tripick,p,t,cc,aa,ff,u,errf,Par);

  if length(i)==0,
    fprintf('\nAdaption completed.\n')
    break;
  elseif size(t,2)>Maxt
    fprintf('\nMaximum number of triangles obtained.\n');
    break
  elseif gen>=Ngen,
    fprintf('\nMaximum number of refinement passes obtained.\n');
    break
  end

  tl=i';
% Kludge: tl must be a column vector
  if size(tl,1)==1,
    tl=[tl;tl];
  end

  u=reshape(u,np,length(u)/np);
  [p,e,t,u]=refinemesh(g,p,e,t,u,tl,Rmethod);
  u=u(:);
  np=size(p,2);
  gen=gen+1;
end


                                                                                                                                                                                                                                                                                       Ned1_0.m                                                                                            0000700 0116504 0001751 00000000634 11757300770 011577  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function shapef=Ned1_0(x,y)
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

return                                                                                                    our_matrices.m                                                                                      0000700 0116504 0001751 00000006741 11757300770 013272  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [A11,M11,L1,S11,Q11,B11,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine assemble and collects several finite element matrices
% given a mesh p,e,t,te. The matrices take into account homogeneous
% boundary condition according to list_of_int_edges and list_of_bnd_edges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% A11 : stiffness matrix for Nedelec elements
% M11 : mass matrix for Nedelec elements
% L11 : load vector with according to function in data.m
% S11 : stiffness matrix for piecewise linear Laplacian
% Q11 : mass matrix for piecewise linear Laplacian
% B11 : coupling matrix for Nedelec/piecewise linear elements
% u_b : boundary condition for electric field
% p_b : boundary condition for potential
% Adjust data for load vector in data.m!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% dimensions
%
int_edge_num=length(list_of_int_edges);
bnd_edge_num=length(list_of_bnd_edges);
int_node_num =length(list_of_int_nodes);
bnd_node_num = length(list_of_bnd_nodes);
edge_num=int_edge_num+bnd_edge_num;
node_num=int_node_num+bnd_node_num;

% permute interior/boundary edges
edge_index=[list_of_int_edges,list_of_bnd_edges];
node_index=[list_of_int_nodes,list_of_bnd_nodes];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble the matrices (not taking into account b.c.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% stiffness and mass matrix for Nedelec elements
[AA,MM]=assemble_maxwell(p,t,te,edge_num);

% stiffness matrix for scalar Laplacian
[SS]=assemble_scalar_laplacian(p,t);

[QQ]=assemble_scalar_mass(p,t);

% Assemble coupling matrix
[BB]=assemble_coupling(p,t,te,edge_num);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set boundary condition for electric field/potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% only zero boundary conditions are possible at the moment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[u_b]=zeros(bnd_edge_num,1);
[p_b]=zeros(bnd_node_num,1);

% loadvector 
[LL]=assemble_load(p,t,te,edge_num);

% permute matrices and load vectors

edge_index=[list_of_int_edges,list_of_bnd_edges];
node_index=[list_of_int_nodes,list_of_bnd_nodes];

AA=AA(edge_index,edge_index);
MM=MM(edge_index,edge_index);
BB=BB(edge_index,node_index);
SS=SS(node_index,node_index);
LL=LL(edge_index);

% define block matrices

A11=AA(1:int_edge_num,1:int_edge_num);
A12=AA(int_edge_num+1:edge_num,1:int_edge_num);
A21=AA(1:int_edge_num,int_edge_num+1:edge_num);
A22=AA(int_edge_num+1:edge_num,int_edge_num+1:edge_num);

L1=LL(1:int_edge_num);

M11=MM(1:int_edge_num,1:int_edge_num);
M12=MM(int_edge_num+1:edge_num,1:int_edge_num);
M21=MM(1:int_edge_num,int_edge_num+1:edge_num);
M22=MM(int_edge_num+1:edge_num,int_edge_num+1:edge_num);

B11=BB(1:int_edge_num,1:int_node_num);
B12=BB(int_edge_num+1:edge_num,1:int_node_num);
B21=BB(int_edge_num+1,int_node_num+1:node_num);
B22=BB(int_edge_num+1:edge_num,int_node_num+1:node_num);

S11=SS(1:int_node_num,1:int_node_num);
S12=SS(int_node_num+1:node_num,1:int_node_num);
S21=SS(1:int_node_num,int_node_num+1:node_num);
S22=SS(int_node_num+1:node_num,int_node_num+1:node_num);

Q11=QQ(1:int_node_num,1:int_node_num);
Q12=QQ(int_node_num+1:node_num,1:int_node_num);
Q21=QQ(1:int_node_num,int_node_num+1:node_num);
Q22=QQ(int_node_num+1:node_num,int_node_num+1:node_num);


return







                               ph_errors.m                                                                                         0000700 0116504 0001751 00000004402 11757300770 012571  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [err0,err1]=ph_errors(p,t,u)

% Computes L2 and H^1 error for scalar FEM solution u
%
%
% u: FEM solution
% t: connectivity list
% p: point list

%%%%%%%%%%%%%%% initialize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x=p(1,:);                             % x-coordinates   
y=p(2,:);                             % y-coordinates

NEL=length(t);                        % number of elements 

[xnod,ynod,wvol]=quad5;               % quadrature points on triangle

[sh]=shap1(xnod,ynod);                % nodal shape functions at quadrature points   
[dxsh,dysh]=dshap1(xnod,ynod);        % and their derivatives

%%%%%%%%%%%%% LOOP OVER ELEMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%

err0=0;
err1=0;

for n=1:NEL                           % LOOP OVER ELEMENTS 

  % elemental data
  
  x0=x(t(1,n));   
  x1=x(t(2,n));
  x2=x(t(3,n));
  y0=y(t(1,n));
  y1=y(t(2,n));
  y2=y(t(3,n));

  B=[x1-x0,x2-x0;y1-y0,y2-y0];        % Jacobian of elemental map      
  trans=[x0;y0];                      % translation vector 
  Binv=inv(B);                        % inverse
  Jdet=det(B);                        % determinant 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ASSEMBLE ERROR CONTRIBUTIONS

  for k=1:length(wvol)                 % LOOP OVER QUADRATURE POINTS

    [pphys]=B*[xnod(k);ynod(k)]+trans; % physical point
    xp=pphys(1);
    yp=pphys(2); 
    dv=Jdet*wvol(k);                   % volume element  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get exact solutions and derivatives
    [dummy1, dummy2, ue, gradue] = exact(xp,yp);
    dxue = gradue(1,1);                        % derivatives      
    dyue = gradue(2,1); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    %%%%%%%%%%%%%% get DG solution and derivative at quad point %%%%

    duh=[0;0];
    uh=0;   
    for i=1:3
      ii=t(i,n);
      uh=uh+u(ii)*sh(i,k);         
      gradi=transpose(Binv)*[dxsh(i,k);dysh(i,k)];
      duh=duh+u(ii)*gradi;
    end
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% add to error %%%%%%%%%%%%%%%%%%%%

    err0=err0+(uh-ue)^2*dv;
    err1=err1+(duh(1)-dxue)^2*dv+(duh(2)-dyue)^2*dv;

  end                                 % LOOP OVER QUADRATURE POINTS

end                                   % LOOP OVER ELEMENTS

err0=sqrt(err0);
err1=sqrt(err1);





















                                                                                                                                                                                                                                                              quad5.m                                                                                             0000700 0116504 0001751 00000001765 11757300770 011616  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [xnod,ynod,w]=quad5

%  
% gives the quadrature points (xnod,ynod) and weights w 
% for the reference triangle.
% 
% This formula is exact for polynomials of degree 5
%


  xnod=zeros(1,7);
  ynod=zeros(1,7);
  w=zeros(1,7);

  third =0.3333333333333333333333333333;
  a=0.0597158717;
  b=0.4701420641;
  c=0.7974269853;
  d=0.1012865073d0;
   
  w(1)=0.225;
  w(2)=0.1323941527;
  w(3)=0.1323941527;
  w(4)=0.1323941527;
  w(5)=0.1259391805;
  w(6)=0.1259391805;
  w(7)=0.1259391805;
  
  lam=zeros(3,7);

  lam(1,1)=third;
  lam(2,1)=third;
  lam(3,1)=third;
 
  lam(1,2)=a;
  lam(2,2)=b;
  lam(3,2)=b;  
  
  lam(1,3)=b;
  lam(2,3)=b;
  lam(3,3)=a;  

  lam(1,4)=b;
  lam(2,4)=a;
  lam(3,4)=b;

  lam(1,5)=c;
  lam(2,5)=d;
  lam(3,5)=d;

  lam(1,6)=d;
  lam(2,6)=d;
  lam(3,6)=c;

  lam(1,6)=d;
  lam(2,6)=d;
  lam(3,6)=c;
  
  lam(1,7)=d;
  lam(2,7)=c;
  lam(3,7)=d;

  xnod=lam(1,:);
  ynod=lam(2,:);  
  w=0.5*w;










           shap1.m                                                                                             0000700 0116504 0001751 00000000347 11757300771 011607  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [sh]=shap1(x,y)
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


                                                                                                                                                                                                                                                                                         sort_edge_dofs.m                                                                                    0000700 0116504 0001751 00000001117 11757300771 013555  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine sorts out the interior and boundary edges 
% of a mesh (p,e,t)
% 
% list_of_int_edges: list of interior edges
% list_of_bnd_edges: list of boundary edges
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

edges=[1:edge_num];

% identify interior/boundary edges

edge_boundary_index=ismember(edges,be);
list_of_bnd_edges=find(edge_boundary_index);
list_of_int_edges=find(1-edge_boundary_index);

return

   
                                                                                                                                                                                                                                                                                                                                                                                                                                                 sort_nodal_dofs.m                                                                                   0000700 0116504 0001751 00000001076 11757300771 013752  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This routine sorts out the interior and boundary degrees nodes 
% of a mesh (p,e,t)
% 
% list_of_int_nodes: list of interior nodes
% list_of_bnd_nodes: list of boundary nodes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nod_num = length(p);
nodes=[1:nod_num];
bnd_nodes=unique(e(1:2,:));
bnd_index=ismember(nodes,bnd_nodes);

list_of_bnd_nodes=find(bnd_index);
list_of_int_nodes=find(1-bnd_index);

return

                                                                                                                                                                                                                                                                                                                                                                                                                                                                  test.m                                                                                              0000700 0116504 0001751 00000013456 11757300771 011557  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                %
%
% Main for mixed low-frequency Maxwell
% 
% we approximate solution (u,p) of
%
% curl curl u + nabla p = f, div u = 0     on the domain [-1,1]x[-1,1]  
% u x n = 0, p = 0  on boundary
%
% by FE-solution (uh, ph) in Ned(0,1) x P_1.
% Ned(0,1) denotes the FE-space of Nedelec edge-elements of first type and
% lowest degree.
% P_1 denotes the space of scalar nodal pw. linear elements.
%
% We create different types of plots for the computed FE-solution.
%
% Anna Schneebeli, June 02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% h = [0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.15];
ref=3;                      % number of global refinements
L2_error_u = zeros(1,ref);
Hcurl_error = zeros(1,ref);
L2_error_p = zeros(1,ref);
H1_error = zeros(1,ref);
l_num = zeros(1,ref);
DoFs = zeros(1,ref);
    
[p,e,t]=initmesh('squareg','hmax',1.5);
for i=1:ref
    
   clear te be AA BB LL A B L uh ph
    
    [p,e,t]=refinemesh('squareg',p,e,t,'regular');
    
    nod_num = length(p);          % number of nodes
  

    [te, be, e_num]=edges(t,e);
    [AA,BB,LL] = assemble_mixed(p,t,te,e_num);
    CC = zeros(nod_num, nod_num);
    [A,B,C,L] = boundary_cond(AA,BB,CC,LL,be,e);
    S=assemble_laplacian(p,t);
[S_bnd] = boundary_cond_laplacian(S,e);
    
  %h2=1/length(t);
  %h=sqrt(h2);
%       disp('Press any key for S_bnd') 
     %vec1(i)=condest(S_bnd);
       gamma1=normest(A)/normest(B)^2;
       gamma2=length(t);
       gamma3=1/length(t);
	gamma4=1 ;
    %vcond1(i)=condest(A+gamma1*B*B');
    vcond2(i)=cond(full(A+gamma2*B*B'))
    %vcond3(i)=condest(A+gamma3*B*B');
    %vcond4(i)=condest(A+gamma4*B*B');
    
    %vmin1(i)=eigs(A+gamma1*B*B',1,'SM');
    %vmin2(i)=eigs(A+gamma2*B*B',1,'SM');
    %vmin3(i)=eigs(A+gamma3*B*B',1,'SM');
    %vmin4(i)=eigs(A+gamma4*B*B',1,'SM');
    %vmax1(i)=eigs(A+gamma1*B*B',1,'LM');
    %vmax2(i)=eigs(A+gamma2*B*B',1,'LM');
    %vmax3(i)=eigs(A+gamma3*B*B',1,'LM');
    %vmax4(i)=eigs(A+gamma4*B*B',1,'LM');
    %if i>1 
	%i%,
	%vcond1(i)/vcond1(i-1),
	%%vcond2(i)/vcond2(i-1),
	%vcond3(i)/vcond3(i-1),
	%vcond4(i)/vcond4(i-1),
	%%vmin1(i)/vmin1(i-1),
	%vmin2(i)/vmin2(i-1),
	%vmin3(i)/vmin3(i-1),
	%vmin4(i)/vmin4(i-1),
	%vmax1(i)/vmax1(i-1),
	%vmax2(i)/vmax2(i-1),
	%vmax3(i)/vmax3(i-1),
	%vmax4(i)/vmin4(i-1),
    
    %end


% 
%  
%   [S_bnd] = boundary_cond_laplacian(S,e)
%    % solve linear system 
%    S = [A B;B' C];
%    F = [L;zeros(nod_num,1)];
%    sol = S\F;
%norm(sol),pause
%    uh = sol(1:e_num,1);
%    ph = sol(e_num+1:length(sol), 1);
%
%    [L2_error_u(i),Hcurl_error(i)] = uh_errors(p,t,te,uh);
%    [L2_error_p(i),H1_error(i)] = ph_errors(p,t,ph);
%
%    % compute number of elements
%    l_num(i) = length(t);
%    
%    % compute number of degrees of freedom in u
%    DoFs(i) = length(uh);
%
%    disp(['mesh ' num2str(i) ': done!'])
end


break;

figure(2)
    if i>1 i,vec1(i)/vec1(i-1),vec2(i)/vec2(i-1),end


% 
%  
%   [S_bnd] = boundary_cond_laplacian(S,e)
%    % solve linear system 
%    S = [A B;B' C];
%    F = [L;zeros(nod_num,1)];
%    sol = S\F;
%norm(sol),pause
%    uh = sol(1:e_num,1);
%    ph = sol(e_num+1:length(sol), 1);
%
%    [L2_error_u(i),Hcurl_error(i)] = uh_errors(p,t,te,uh);
%    [L2_error_p(i),H1_error(i)] = ph_errors(p,t,ph);
%
%    % compute number of elements
%    l_num(i) = length(t);
%    
%    % compute number of degrees of freedom in u
%    DoFs(i) = length(uh);
%
%    disp(['mesh ' num2str(i) ': done!'])
end


break;

figure(2)
    if i>1 i,vec1(i)/vec1(i-1),vec2(i)/vec2(i-1),end


% 
%  
%   [S_bnd] = boundary_cond_laplacian(S,e)
%    % solve linear system 
%    S = [A B;B' C];
%    F = [L;zeros(nod_num,1)];
%    sol = S\F;
%norm(sol),pause
%    uh = sol(1:e_num,1);
%    ph = sol(e_num+1:length(sol), 1);
%
%    [L2_error_u(i),Hcurl_error(i)] = uh_errors(p,t,te,uh);
%    [L2_error_p(i),H1_error(i)] = ph_errors(p,t,ph);
%
%    % compute number of elements
%    l_num(i) = length(t);
%    
%    % compute number of degrees of freedom in u
%    DoFs(i) = length(uh);
%
%    disp(['mesh ' num2str(i) ': done!'])
end


break;

figure(2)
subplot(2,1,1)
loglog(l_num,L2_error_u,'o-',l_num,Hcurl_error,'*-',l_num,l_num.^(-1/2),'r')
% loglog(l_num,Hcurl_error,'*-',l_num,l_num.^(-1/2),'r')
title('L2- and H(curl)-error of uh vs. number of elements')
legend('L2-error','H(curl)-error','Order 1')
%title('H(curl)-error vs. number of elements')
%legend('H(curl)-error','Order 1')
xlabel('number of elements in mesh')
ylabel('norm-errors')
grid on

subplot(2,1,2)
loglog(l_num,L2_error_p,'o-',l_num,H1_error,'*-',l_num,l_num.^(-1/2),l_num,l_num.^(-1),'-.' )
% loglog(l_num,Hcurl_error,'*-',l_num,l_num.^(-1/2),'r')
title('L2- and H1-error of ph vs. number of elements')
legend('L2-error','H1-error','Order 1', 'Order 2')
%title('H(curl)-error vs. number of elements')
%legend('H(curl)-error','Order 1')
xlabel('number of elements in mesh')
ylabel('norm-errors')
grid on



% plots of the solution (uh, ph)
%
x=linspace(-1,1,15);
y=x;

% surface plot of uh
u_fem = uh_on_grid(p,t,te,uh,x,y);
disp('here we go! 2-D plots...')
z1=zeros(length(y),length(x));
z2=z1;
z1(:,:)=u_fem(1,:,:);
z2(:,:)=u_fem(2,:,:);

% plot of x-component of uh
figure(3)
% pdemesh(p,e,t,zeros(1,length(p)))
% hold on
mesh(x,y,z1)
title('x-component of FE-solution')
xlabel('x')
ylabel('y')
hidden on

% plot of y-component of uh
figure(4)
surf(x,y,z2)
title('y-component of FE-solution')
xlabel('x')
ylabel('y')


% vector-field plot of uh
disp('here we go! vector plots...')
figure(5)
quiver(x,y,z1,z2)
title('Vector-field plot of FE-solution uh')
xlabel('x')
ylabel('y')

% plot of the solution ph on finest grid with pde-toolbox commands
figure(6)
pdeplot(p,e,t, 'xydata', ph, 'zdata', ph, 'mesh', 'on')
title('FE-pressure ph')


% Here we automatically save the data in a file called [filename].mat
% The data can be extracted by typing the command >> load [filename]  in the MATLAB command window
%save example2

%return
                                                                                                                                                                                                                  uh_errors.m                                                                                         0000700 0116504 0001751 00000004575 11757300771 012612  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function [L2, Hcurl] = uh_errors(p,t,te,uh)
%
% Anna Schneebeli, June 02
%
% compute L^2- and H(curl)-norm of the error u-uh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('compute L2- and H(curl)-error....')

L2 = 0;
Hcurl = 0;

% get point and element data 
%
el_num =length(t);                          % number of elements
x=p(1,:);                                  % global x-coordinates of vertices
y=p(2,:);                                  % global y-coordinates of vertices

% quadrature points on master-element
%
[xnod,ynod,weight]=quad5;

% get reference shape functions[xnod,ynod,weight]=quad5;
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
 
  B=[x1-x0,x2-x0;y1-y0,y2-y0];        % Jacobian of affine element-map      
  trans=[x0;y0];                      % translation vector of affine element-map 
  Binv=inv(B);                        % inverse
  Jdet=det(B);                        % determinant of B
  BinvT=transpose(Binv);              % transposed of Binv
  dvol = Jdet * weight;               % volume element
  
  % get edge-dof numbers of current element
  [dof] = te(1:3,l);
        
  % get orientation of edges
  [orient] = diag(te(4:6,l));
        
  % get uh on dofs of current element, orientation inclusive
  [uhl] = orient*uh(dof);
  
  for k=1:length(weight)        % loop over quadrature points
      
    [quad]=B*[xnod(k);ynod(k)]+trans;                         % quadrature points mapped onto actual element
    [u, curlu] = exact(quad(1),quad(2));            % exact solution and curl of it on transformed quadrature points
   
    % u_fem on current element at quadrature point (Piola-Trafo of u_fem on master element at quad. pts.)
    u_fem = BinvT * sh_ref(:,:,k)*uhl;
    
    % curl u_fem on current element at quadrature point; curl(sh_ref)=2, and curl(sh)=1/Jdet*curl(sh_ref) for every shapefkt.
    curlu_fem = 2/Jdet*sum(uhl);
    
    % compute L2-error (square)
    L2 = L2 + norm((u(:,1)-u_fem))^2*dvol(k);
    
    % compute Hcurl-seminorm (square) of error 
    Hcurl = Hcurl + (curlu-curlu_fem)^2*dvol(k);
    
  end           % loop over quadrature point
end             % loop over elements

L2 = sqrt(L2);
Hcurl = sqrt(Hcurl);

return
  
  
  


                                                                                                                                   uh_on_grid.m                                                                                        0000700 0116504 0001751 00000007120 11757300771 012704  0                                                                                                    ustar   greif                           faculty                                                                                                                                                                                                                function u_fem = uh_on_grid(P,T,Te,uh,x,y)
%
% Anna Schneebeli, June 02
%
% computes the Nedelec(1,0)-FE-solution u_fem at gridpoints (x(i),y(j)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Nedelec(1,0) means that we use Nedelec shape-functions of first type and lowest order to approximate u)
%
% This function enables to visualize the computed FE-solution componentwise on the
% plot-grid [x,y] using the commands 
% 
% Z=zeros(length(x),length(y));
% Z(:,:)= u_fem(1,:,:);      respectively:  Z(:,:)=,u_fem(2,:,:);
% surf(x,y,Z);
%
% TODO: make sure, that dofs on dirichlet-boundary have the correct value!! 
% This is not guaranteed on points which coincide with boundary-vertices that belong to more than 
% one element. This stems in the fact that the tria2grid routine is made for nodal FEM and there, we have not
% to care about which of the neighbouring element a vertex belongs to!
% Conversely, in edge-FEM, the nodal value depends on the edge from which a vertex is looked at!
% 
% Remark: In general, the x- and y-components of the FE-solution are not continuous over edges! Tri2grid assigns 
% points to a unique element via "forward assignment" (i.e. to the element right to the edge).

disp('get data for plotting...')

nx = length(x);
ny = length(y);
u_fem = zeros(2,ny,nx);

% get element indices of grid points
udummy = zeros(length(P(1,:)),1);               % uh are values at edge-DOFs, but PDEtoolbox-function tri2grid needs nodal dofs,
                                                % so here they are created artificially....
[dummy, ptind] = tri2grid(P,T,udummy,x,y);

% loop over all grid points
for j = 1:ny
    for i = 1:nx
        
        % get element number l of point (x(i),y(j)) 
        % !! CAREFUL: the x-values correspond to the COLUMNS of ptind and thus the index order of i and j is inverted !!  
        l = ptind(j,i);
        
        % get global coordinates of vertices 1, 2, 3 of current element
        x0=P(1,(T(1,l)));   
        x1=P(1,(T(2,l)));
        x2=P(1,(T(3,l)));
        y0=P(2,(T(1,l)));
        y1=P(2,(T(2,l)));
        y2=P(2,(T(3,l)));
       
        B=[x1-x0,x2-x0;y1-y0,y2-y0];        % Jacobian of affine element-map
        trans=[x0;y0];                      % translation vector
        Binv = inv(B);
        
       
        % get values of (transformed) shape functions at grid-point (x(i),y(j))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pull back point
        refpt = Binv*([x(i);y(j)]-trans);

        % reference shape functions
        shapef = Ned1_0(refpt(1),refpt(2));
        
        % push forward shape functions (Piola-Trafo)
        lshapef = transpose(Binv)*shapef(:,:,1);
        
        % get edge-dof numbers of current element
        dof = Te(1:3,l);
        
        % get orientation of edges
        orient = diag(Te(4:6,l));
        
        
        % assemble u_fem at grid point (x(i),y(j))
        %
        % !! CAREFUL: the x-values correspond to the COLUMNS of u_fem(comp,:,:) and thus the index order of i and j is inverted !! 
        % !! (The MATLAB plot-functions 'surf' and 'mesh' require this definition of the rows and columns of u_fem) !!
        %
        u_fem(:,j,i) = u_fem(:,j,i)+ lshapef*orient*uh(dof);
 
%         helpu = exact(x(i),y(j));
%         u_fem(1,j,i) = helpu(1);
%         u_fem(2,j,i) = helpu(2);
        
         % correction for boundary-vertices (1,-1),(-1,1) and (1,1)
        if (x(i)==1 & y(j)==-1) | (x(i)==1 & y(j)==1) | (x(i)==-1 & y(j)==1)
            u_fem(1,j,i)=0;
            u_fem(2,j,i)=0;
        end
      
    end
end

return
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                