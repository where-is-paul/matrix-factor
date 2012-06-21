//-*-mode:c++-*-
#ifndef _TRIPLET_MATRIX_GAXPY_H_
#define _TRIPLET_MATRIX_GAXPY_H_

template<class idx_type, class el_type>
void triplet_matrix<idx_type, el_type> :: gaxpy (const el_type* x, el_type* y) const
{
    idx_type nnz_ = nnz ();
    for (idx_type elem = 0; elem < nnz_; elem ++)
    {
        y[m_row_idx[elem]] += m_x[elem] * x[m_col_idx[elem]];
    }
}

template<class idx_type, class el_type>
void triplet_matrix<idx_type, el_type> :: gaxpy (const std::vector<el_type>& x, std::vector<el_type>& y) const
{
    if (!(n_cols () == y.size() && n_rows() == x.size()))
    {
        throw std::logic_error(__FUNCTION__);
    }

    gaxpy (&(x[0]), &(y[0]));
}
#endif // _TRIPLET_MATRIX_GAXPY_H_
