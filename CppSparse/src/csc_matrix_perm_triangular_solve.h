// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_PERM_TRIANGULAR_SOLVE_H_
#define _CSC_MATRIX_PERM_TRIANGULAR_SOLVE_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: lperm (idx_vector_type& perm) const
{
    perm.resize(m_n_cols);
    std::fill(perm.begin(), perm.end(), npos);

    for (idx_type col = m_n_cols - 1; col >= 0; col --)
    {
        idx_type num_marked_rows = 0;
        idx_type marked_offset   = npos;

        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            if (perm[row] != npos)
            {
                num_marked_rows ++;
                marked_offset = offset;
            }            
        }
        if (num_marked_rows != 1)
        {
            throw std::logic_error("Matrix is not permuted lower-triangular");
        }
        perm[col] = marked_offset;
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: uperm (idx_vector_type& perm) const
{
    perm.resize(m_n_cols);
    std::fill(perm.begin(), perm.end(), npos);

    for (idx_type col = 0; col < m_n_cols; col ++)
    {
        idx_type num_marked_rows = 0;
        idx_type marked_offset   = npos;

        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            if (perm[row] != npos)
            {
                num_marked_rows ++;
                marked_offset = offset;
            }
        }

        if (num_marked_rows != 1)
        {
            throw std::logic_error("Matrix is not permuted upper-triangular");
        }

        perm[col] = marked_offset;
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: perm_lsolve  (elt_vector_type& b) const
{
    idx_vector_type perm(m_n_cols);
    lperm(perm);
    
    /**
     *  Back-solve:
     *  for j = 0 : n-1
     *     b[j] /= A[j, j];
     *     for i = j + 1 : n - 1
     *         b[i] -= A[i, j] * b[j]
     *     end
     *  end
     */
    for (idx_type j = 0; j < m_n_cols; j ++)
    {
        idx_type diag_offset = perm[j];
        el_type& bj          = b[j];

        // First compute x[j]
        bj /= m_x[diag_offset];
        assert(diag_offset < m_col_idx[j + 1]);

        idx_type offset = m_col_idx[j];
        for (; offset < diag_offset; offset ++)
        {
            b[m_row_idx[offset]] -= m_x[offset] * bj;
        }
        // jump over the diagonal element
        offset ++;
        // continue to the rest of the rows
        for (; offset < m_col_idx[j + 1]; offset ++)
        {
            b[m_row_idx[offset]] -= m_x[offset] * bj;
        }
    }
}

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: perm_usolve  (elt_vector_type& b) const
{
    idx_vector_type perm(m_n_cols);
    uperm(perm);
    
    /**
     * Back-solve:
     * for j = n-1 : 0
     *   b[j] /= A[j, j]
     *   for i = 0 : j - 1
     *      b[i] -= A[i, j] * b[j]
     *   end
     * end
     */

    for (idx_type j = m_n_cols - 1; j >= 0; j --)
    {
        idx_type diag_offset = perm[j];
        el_type& bj          = b[j];
        b[j] /= m_x[diag_offset];

        idx_type offset = m_col_idx[j];
        for (; offset < diag_offset; offset ++)
        {
            b[m_row_idx[offset]] -= m_x[offset] * bj;
        }
        offset ++;
        for (; offset < m_col_idx[j + 1]; offset ++)
        {
            b[m_row_idx[offset]] -= m_x[offset] * bj;
        }
        
    }
}


#endif
