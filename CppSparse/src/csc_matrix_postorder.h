// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_POSTORDER_
#define _CSC_MATRIX_POSTORDER_

template<class idx_type, class el_type>
typename csc_matrix<idx_type, el_type>::idx_vector_type 
csc_matrix<idx_type, el_type>::postorder(const idx_vector_type& parent) const
{
    // The first-child, next-sibling terminology is based on 
    // how ANTLR names its trees. 
    idx_vector_type first_child(m_n_rows, npos), next_sibling(m_n_rows, npos);

    // Inspect the parent of each node
    // If node is a root, skip it
    // If not, add it to the previous sibling of the parent's last
    // child and make it the first child of the parent. This ensures that
    // the first child of a node is the one with the highest number, i.e.,
    // the list of children are sorted in reverse order
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        if (parent[col] != npos)
        {
            auto sibling             = first_child[parent[col]];
            first_child[parent[col]] = col;
            next_sibling[col]        = sibling;            
        }
    }

    idx_vector_type stack;
    idx_vector_type post;
    
    stack.reserve(m_n_cols); // Pushs and Pops can be constant time as opposed to amortized constant time.
    post.reserve(m_n_cols);

    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        if (parent[col] == npos)
        {
            dfs_tree (col, first_child, next_sibling, post, stack);
        }
    }

    return post;
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type>::dfs_tree (idx_type col, 
                                              idx_vector_type& first_child, 
                                              const idx_vector_type& next_sibling, 
                                              idx_vector_type& post, 
                                              idx_vector_type& stack) const
{
    stack.push_back(col);
    while(!stack.empty())
    {
        idx_type vertex             = stack.back();
        idx_type first_child_vertex = first_child[vertex];
        // If this vertex has no children then just push it to the
        // postordered tree
        if (first_child_vertex == npos)
        {
            stack.pop_back();
            post.push_back(vertex);
        }
        else
        {
            stack.push_back(first_child_vertex);
            first_child[vertex] = next_sibling[first_child_vertex];
        }
    }
}


#endif
