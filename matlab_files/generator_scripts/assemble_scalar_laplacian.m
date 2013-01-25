function [A]=assemble_scalar_laplacian(p,t)

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



































