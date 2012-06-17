// -*- mode: c++ -*-

#ifndef _CSC_MATRIX_SP_TRIANGULAR_SOLVE_H_
#define _CSC_MATRIX_SP_TRIANGULAR_SOLVE_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type>::lsolve (const csc_t& b, idx_vector_type& x, idx_vector_type& xi, idx_type k /* = 0 */) const
{
    if (xi.size() < n_rows())
    {
        throw std::runtime_error("Insuffient length in work vector.");
    }

    idx_vector_type pinv;
    idx_type top = reach(b, k, pinv, xi);

    for (idx_type offset = top; offset < m_n_cols; offset ++)
    {
        x[xi[offset]] = 0;
    }

    for (idx_type offset = b.m_col_idx[k]; offset < b.m_col_idx[k + 1]; offset ++)
    {
        x[b.m_row_idx[offset]] = b.m_x[offset];
    }

    auto col_begin = xi.begin() + top;
    auto col_end   = xi.begin() + n_cols();
    for (auto col_iter = col_begin; col_iter != col_end; col_iter ++)
    {
        idx_type col = *col_iter;
        idx_type offset = m_col_idx[col];
        assert(col == m_row_idx[offset]);
        x[col] /= m_x[offset];
        offset ++;
        for (; offset < m_col_idx[col + 1]; offset ++)
        {
            x[m_row_idx[offset]] -= x[col] * m_x[offset];
        }
    }
}

#endif
