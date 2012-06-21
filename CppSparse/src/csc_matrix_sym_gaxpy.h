//-*-mode:c++-*-
#ifndef _CSC_MATRIX_SYM_GAXPY_H_
#define _CSC_MATRIX_SYM_GAXPY_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: sym_gaxpy (const el_type* x, el_type* y, uplo_type type) const
{
    if (type == UPPER_TRIANGULAR)
    {
        for (idx_type col = 0; col < n_cols(); col ++)
        {
            for (idx_type offset = m_col_idx[col]; offset < m_col_idx [col + 1]; offset ++)
            {
                idx_type row = m_row_idx [offset];
                el_type  elt = m_x       [offset];
                if (row < col)
                {
                    y[row] += elt * x[col];
                    y[col] += elt * x[row];
                }
                else if (row == col)
                {
                    y[row] += elt * x[col];
                }
                else
                {
                    break;
                  }
            }
        }
    }
    else
    {
        for (idx_type col = 0; col < n_cols(); col ++)
        {
            for (idx_type offset = m_col_idx[col]; offset < m_col_idx [col + 1]; offset ++)
            {
                idx_type row = m_row_idx [offset];
                if (row < col)
                {
                    continue;
                }
                else if (row == col)
                {
                    el_type  elt = m_x [offset];
                    y[row] += elt * x[col];
                }
                else
                {
                    el_type  elt = m_x [offset];
                    y[row] += elt * x[col];
                    y[col] += elt * x[row];
                }
            }
        }
    }
}

template<class idx_type, class el_type>
void 
csc_matrix<idx_type, el_type> :: sym_gaxpy (const std::vector<el_type>& x, std::vector<el_type>& y, uplo_type side) const
{
    if (!(n_rows() == n_cols() && n_rows() == x.size() && n_rows () == y.size()))
    {
        throw std::logic_error(__FUNCTION__);
    }

    sym_gaxpy (&(x[0]), &(y[0]), side);
}

#endif // _CSC_MATRIX_SYM_GAXPY_H_
