function u_fem = uh_on_grid(P,T,Te,uh,x,y)
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
