//-*-mode:c++-*-
#ifndef _CSC_MATRIX_GAXPY_H_
#define _CSC_MATRIX_GAXPY_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: gaxpy (const el_type* x, el_type* y) const
{
    for (idx_type col = 0; col < n_cols(); col ++)
    {
        for (idx_type row = m_col_idx[col]; row < m_col_idx[col + 1]; row ++)
        {
            y [m_row_idx[row]] += m_x [row] * x[col];
        }
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: gaxpy (const std::vector<el_type>& x, std::vector<el_type>& y) const
{
    if (!(n_cols () == y.size() && n_rows() == x.size()))
    {
        throw std::logic_error(__FUNCTION__);
    }

    gaxpy (&(x[0]), &(y[0]));
}

#endif // _CSC_MATRIX_GAXPY_H_
