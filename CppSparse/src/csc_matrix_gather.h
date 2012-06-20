// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_GATHER_H_
#define _CSC_MATRIX_GATHER_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: gather (idx_type column, const elt_vector_type& vec)
{
    for (idx_type offset = m_col_idx[column]; offset < m_col_idx[column+1]; offset ++)
    {
        m_x[offset] = vec[m_row_idx[offset]];
    }
}

#endif
