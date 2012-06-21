// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_PERM_TRIANGULAR_SOLVE_GENERAL_H_
#define _CSC_MATRIX_PERM_TRIANGULAR_SOLVE_GENERAL_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: triangular_solve (elt_vector_type& b) const
{
    if (m_n_rows != m_n_cols)
    {
        throw std::logic_error("Matrix must be square");
    }

    idx_type N = m_n_rows;

    // row_counts corresponds to the "r" array in the book
    idx_vector_type row_counts(N);
    // col_indices corresponds to the "z" array in the book
    idx_vector_type col_indices(N);

    for (idx_type col = 0; col < N; col ++)
    {
        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type row = m_row_idx[offset];
            row_counts [row] ++;
            col_indices[row] ^= col;
        }
    }
    
    idx_vector_type singletons;
    singletons.reserve(N);
    for (idx_type row = 0; row < N; row ++)
    {
        if (row_counts[row] == 1)
        {
            singletons.push_back(row);
        }
    }

    idx_vector_type perm_rows(N);
    idx_vector_type perm_cols(N);
    idx_vector_type perm_offs(N);
    
    for (idx_type iter = 0; iter < N; iter ++)
    {
        if (singletons.empty())
        {
            throw std::logic_error("Matrix is not permuted upper or lower-triangular");
        }

        idx_type row = singletons.back();
        idx_type col = col_indices[row];

        singletons.pop_back();
        
        perm_rows[iter] = row;
        perm_cols[iter] = col;

        for (idx_type offset = m_col_idx[col]; offset < m_col_idx[col + 1]; offset ++)
        {
            idx_type t_row = m_row_idx[offset];
            col_indices [t_row] ^= col;
            row_counts [t_row] --;
            if (row_counts[t_row] == 1)
            {
                singletons.push_back(t_row);
            }
            if (t_row == row)
            {
                perm_offs[iter] = offset;
            }
        }
    }
    std::cout << perm_rows << std::endl;
    std::cout << perm_cols << std::endl;
    std::cout << perm_offs << std::endl;

    for (idx_type pcol = 0; pcol < m_n_cols; pcol ++)
    {
        idx_type j           = perm_cols[pcol];
        el_type& bj          = b[j];
        idx_type diag_offset = perm_offs[j];
        el_type  ajj         = m_x[diag_offset];
        bj                  /= ajj;

        idx_type offset;
        for (offset = m_col_idx[j]; offset < diag_offset; offset ++)
        {
            idx_type i    = perm_offs[m_row_idx[offset]];
            el_type& bi   = b[i];
            el_type  aij  = m_x[offset];
            bi           -= aij * bj;
        }
        // Skip over the diagonal row
        offset ++;

        for (offset; offset < m_col_idx[j + 1]; offset ++)
        {
            idx_type i    = perm_offs[m_row_idx[offset]];
            el_type& bi   = b[i];
            el_type  aij  = m_x[offset];
            bi           -= aij * bj;
        }
    }

}


#endif
