//-*-mode: c++-*-
#ifndef _CSC_MATRIX_ETREE_H_
#define _CSC_MATRIX_ETREE_H_

template<class idx_type, class el_type>
typename csc_matrix<idx_type, el_type>::idx_vector_type csc_matrix<idx_type, el_type>::etree(bool trans) const
{
    if (n_rows() != n_cols() && !trans)
    {
        throw std::logic_error("Matrix must be square.");
    }

    auto m = n_rows(), n = n_cols();
    // parent[I] is the immediate parent of node I 
    // in the elimination tree.
    idx_vector_type parent(n, npos);
    idx_vector_type ancestor(n, npos);
    idx_vector_type prev(m, npos);

    for (idx_type col = 0; col < n; col ++)
    {
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row  = m_row_idx[offset];
            idx_type i    = trans ? prev[row] : row;
            idx_type next_ancestor = npos;
            for (; i != npos && i < col; i = next_ancestor)
            {
                next_ancestor = ancestor[i];
                // Highest ancestor of i is always col since L[k,i] != 0 <=> A[i, k] != 0
                ancestor[i]   = col;
                // If i is a root node of its sub-tree
                // then col is its parent
                // -- this is referred to as "path compression"
                if (parent[i] == npos)
                {
                    parent[i] = col;
                }
            }
            if (trans)
            {
                prev[row] = col;
            }
        }
    }
    return parent;
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type>::etree_plot (std::ostream& stream, bool trans) const
{
    std::string label = "A";
    if (trans)
    {
        label = "ATA";
    }
        
    stream << "digraph " << label << "{ " << std::endl;
    auto parent = etree(trans);

    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        if (parent[col] != npos)
        {
            // NOTE: The way the DOT language is specified, to indicate a tree you always
            // NOTE: draw the directed edge from PARENT -> CHILD as opposed to
            // NOTE: CHILD -> PARENT which is the convention used to represent the
            // NOTE: elimination tree.
            stream << parent[col] << " -> " << col << " [dir = none];" << std::endl;
            // This creates an "upside-down" elimination tree
            // stream << col << " -> " << parent[col] << " [dir = none]; " << std::endl;
        }
    }


    stream << "}" << std::endl;
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type>::etree_plot (const std::string& file_name, bool trans) const
{
    std::ofstream output_stream(file_name.c_str());
    etree_plot(output_stream);
}

template<class idx_type, class el_type>
std::string csc_matrix<idx_type, el_type>::etree_plot (bool trans) const
{
    std::ostringstream os;
    etree_plot(os);
    return os.str();
}

#endif
