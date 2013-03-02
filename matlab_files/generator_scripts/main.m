%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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









