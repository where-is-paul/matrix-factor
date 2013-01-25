function [err0,err1]=ph_errors(p,t,u)

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





















