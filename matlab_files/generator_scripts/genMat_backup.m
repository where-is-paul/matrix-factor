%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN PROGRAM TO COMPUTE OUR MATRICES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:1,
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








