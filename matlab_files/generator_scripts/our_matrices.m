function [A11,M11,L1,S11,Q11,B11,u_b,p_b]=our_matrices(p,t,te,list_of_int_edges,list_of_bnd_edges,list_of_int_nodes,list_of_bnd_nodes)

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







