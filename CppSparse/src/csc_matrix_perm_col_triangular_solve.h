// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_PERM_COL_TRIANGULAR_SOLVE_H_
#define _CSC_MATRIX_PERM_COL_TRIANGULAR_SOLVE_H_

#include <algorithm>

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: lower_col_perm (idx_vector_type& perm) const
{
    std::fill(perm.begin(), perm.end(), npos);
    for (idx_type col = 0; col < m_n_cols - 1; col ++)
    {
        idx_type row = m_row_idx[m_col_idx[col]];
        perm[row] = col;
    }
    if (std::any_of(perm.cbegin(), perm.cend(), [](idx_type perm){return perm == npos;}))
    {
        throw std::logic_error("Matrix is not column permuted lower-triangular.");
    }
}


template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: upper_col_perm (idx_vector_type& perm) const
{
    std::fill(perm.begin(), perm.end(), npos);

    for (idx_type col = 0; col < m_n_cols - 1; col ++)
    {
        // The diagonal element of each column is the last
        // row in that column
        idx_type offset = m_col_idx[col + 1] - 1;
        idx_type row    = m_row_idx[offset];
        perm[row]       = col;
    }

    if (std::any_of(perm.cbegin(), perm.cend(), [](idx_type perm){return perm == npos;}))
    {
        throw std::logic_error("Matrix is not column permuted upper-triangular.");
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: lower_col_perm_solve  (elt_vector_type& b) const
{
    idx_vector_type perm(m_n_cols);
    lower_col_perm(perm);        

    for (idx_type pcol = 0; pcol < m_n_cols; pcol ++)
    {
        idx_type col = perm[pcol];
        el_type& bj  = b[col];

        // For the diagonal element, compute b[j] /= A[j,j]
        bj /= m_x[m_col_idx[col]];
        
        // For the off-diagonal elements, compute b[i] -= A[i, j] x b[j]
        for (idx_type offset = m_col_idx[col] + 1; offset < m_col_idx[col + 1]; offset ++)
        {
            b[m_row_idx[offset]] -= m_x[offset] * bj;
        }
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: upper_col_perm_solve  (elt_vector_type& b) const
{
    idx_vector_type perm (m_n_cols);
    upper_col_perm (perm);

    for (idx_type pcol = m_n_cols - 1; pcol >= 0; pcol --)
    {
        idx_type col         = perm[pcol];
        idx_type diag_offset = m_col_idx[col+1] - 1;
        el_type& bj          = b[col];

        // b[j] = b[j] / A[j, j]
        bj /= m_x[diag_offset];

        for (idx_type offset = m_col_idx[col]; offset < diag_offset; offset ++)
        {
            // b[i] -= A[i, j] * b[j] 
            b[m_row_idx[offset]] -= m_x[offset] * bj;
        }
    }
}

#endif
