function [L2, Hcurl] = uh_errors(p,t,te,uh)
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
  
  
  


