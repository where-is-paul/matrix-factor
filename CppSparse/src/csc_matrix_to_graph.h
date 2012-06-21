// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_TO_GRAPH_H_
#define _CSC_MATRIX_TO_GRAPH_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: to_graph (std::ostream& stream, bool trans) const
{
    stream << "digraph " << (trans ? std::string("G_T") : std::string("G")) << " { " << std::endl ;
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {           
            idx_type i = m_row_idx[offset], j = col;
            if (i != j)
            {
                if (trans)
                {
                    stream << j << " -> " << i << "; " << std::endl;
                }
                else
                {
                    stream << i << " -> " <<  j << "; " << std::endl;
                }
            }
        }
    }
    stream << "}" << std::endl;
}


template<class idx_type, class el_type>
std::string csc_matrix<idx_type, el_type> :: to_graph (bool trans) const
{
    std::ostringstream os;
    to_graph(os, trans);
    return os.str();
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: to_graph (const std::string& file_name, bool trans) const
{
    std::ofstream file(file_name.c_str());
    to_graph(file, trans);
}

#endif // _CSC_MATRIX_TO_GRAPH_H_
