//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_FILTER_H_
#define _CSC_MATRIX_FILTER_H_

template<class idx_type, class el_type>
template<class predicate>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: filter (const predicate& pred) const
{
    csc_t ret(m_n_rows, m_n_cols, nnz());
    
    idx_type current_nnz = 0;
    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        ret.m_col_idx[col] = current_nnz;
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx [offset];
            el_type  val = m_x       [offset];
            if (pred(row, col, val))
            {
                ret.m_row_idx.push_back(row);
                ret.m_x.push_back(val);
                current_nnz ++;
            }
        }
    }
    ret.m_col_idx.back() = current_nnz;
    ret.m_row_idx.resize(current_nnz);
    ret.m_x.resize(current_nnz);
    return ret;    
}

#endif
