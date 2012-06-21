//-*-mode:c++-*-
#ifndef _CSC_MATRIX_OPERATOR_EQ_H_
#define _CSC_MATRIX_OPERATOR_EQ_H_

template<class idx_type, class el_type>
bool csc_matrix<idx_type, el_type> :: operator== (const csc_matrix<idx_type, el_type>& other) const
{
    // Compare the sizes and number of non-zero entries
    if (!((n_rows() == other.n_rows()) &&
          (n_cols() == other.n_cols()) &&
          (nnz()    == other.nnz ())))
    {
        return false;
    }

    if (!(m_col_idx == other.m_col_idx &&
          m_row_idx == other.m_row_idx &&
          m_x       == other.m_x))
    {
        return false;
    }

    return true;
}

#endif
