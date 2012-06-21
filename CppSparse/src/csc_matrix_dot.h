//-*-mode:c++-*-
#ifndef _CSC_MATRIX_DOT_H_
#define _CSC_MATRIX_DOT_H_

template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: dot_is_sorted (const csc_matrix<idx_type, el_type>& y) const
{
    const csc_matrix<idx_type, el_type>& x = *this;

    idx_type offset_x = 0;
    idx_type offset_y = 0;
    idx_type nnz_x    = x.nnz();
    idx_type nnz_y    = y.nnz();
    el_type  result   = static_cast<el_type>(0);

    while (offset_x < nnz_x && offset_y < nnz_y)
    {
        if (x.m_row_idx[offset_x] == y.m_row_idx[offset_y])
        {
            result += x.m_x[offset_x] * y.m_x[offset_y];
            offset_x ++;
            offset_y ++;
        }
        else if (x.m_row_idx[offset_x] < y.m_row_idx[offset_y])
        {
            offset_x ++;
        }
        else
        {
            offset_y ++;
        }
    }
    return result;
}
    

template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: dot_no_sorted (const csc_matrix<idx_type, el_type>& y) const
{
    const csc_matrix<idx_type, el_type>& x = *this;
    
    idx_type nnz_x    = x.nnz();
    idx_type nnz_y    = y.nnz();
    el_type  result   = static_cast<el_type>(0);
    
    std::vector<el_type> work (x.n_rows());
    
    // Scatter X into the work vector
    for (idx_type offset = 0; offset < nnz_x; offset ++)
    {
        work[x.m_row_idx[offset]] = x.m_x[offset];
    }

    // Compute dot product of work and Y    
    for (idx_type offset = 0; offset < nnz_y; offset ++)
    {
        result += work[y.m_row_idx[offset]] * y.m_x[offset];
    }
                               
                               

    return result;
}


template<class idx_type, class el_type>
el_type csc_matrix<idx_type, el_type> :: dot (const csc_matrix<idx_type, el_type>& y, bool is_sorted) const
{
    if (!(n_cols() == 1 && y.n_cols() == 1 && n_rows() == y.n_rows()))
    {
        throw std::logic_error(__FUNCTION__);
    }

    if (is_sorted)
    {
        return dot_is_sorted (y);
    }
    else
    {
        return dot_no_sorted (y);
    }
}




#endif // _CSC_MATRIX_DOT_H_
