// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_PERM_H_
#define _CSC_MATRIX_PERM_H_

template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: col_permute (const idx_vector_type& q) const
{
    csc_t c (n_rows(), n_cols(), nnz());
    
    idx_type nnz = 0, offset = 0;
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        idx_type a_col = q[col];
        for (idx_type offset = m_col_idx[a_col]; offset < m_col_idx[a_col + 1]; offset ++)
        {
            c.m_row_idx[nnz] = m_row_idx[offset];
            c.m_x      [nnz] = m_x      [offset];
            nnz ++;
        }
        c.m_col_idx[col+1] = nnz;
    }
    return c;
}


template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: row_permute (const idx_vector_type& p) const
{
    csc_t c (n_rows(), n_cols(), nnz());
    
    std::copy (m_col_idx.begin(), m_col_idx.end(), c.m_col_idx.begin());
    for (idx_type offset = 0; offset < nnz(); offset ++)
    {
        c.m_row_idx[offset] = m_row_idx[p[offset]];
    }
    return c;
}



template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: permute (const idx_vector_type& p, perm_type type) const
{
    switch (type)
    {
    case PERM_ROW:
        return row_permute(p);
    case PERM_COL:
        return col_permute(p);
    case PERM_SYM:
        return sym_permute(p);
    default:
        throw std::logic_error("Invalid permutation type.");
    }
}

template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: permute (const idx_vector_type& p, const idx_vector_type& q) const
{
    csc_t c(n_rows(), n_cols(), nnz());
    idx_type nnz = 0;
    for (idx_type col = 0; col < n_cols(); col ++)
    {
        idx_type pcol = q[col];
        for (idx_type offset = m_col_idx[pcol]; offset < m_col_idx[pcol+1]; offset ++)
        {
            c.m_row_idx[nnz] = p[m_row_idx[offset]];
            c.m_x      [nnz] = m_x[nnz];
            nnz ++;
        }
        c.m_col_idx[col + 1] = nnz;
    }
    return c;
}

#endif
