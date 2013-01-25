function [list_of_int_edges,list_of_bnd_edges]=sort_edge_dofs(te,be,edge_num)
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

   
