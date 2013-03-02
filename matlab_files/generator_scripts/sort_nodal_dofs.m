function [list_of_int_nodes,list_of_bnd_nodes]=sort_nodal_dofs(p,e)
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

