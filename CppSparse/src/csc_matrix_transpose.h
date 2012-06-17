//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_TRANSPOSE_H_
#define _CSC_MATRIX_TRANSPOSE_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: transpose(csc_matrix<idx_type, el_type>& result) const
{
    result.m_n_rows = m_n_cols;
    result.m_n_cols = m_n_rows;

    const idx_vector_type& row_idx = m_row_idx;
    const idx_vector_type& col_idx = m_col_idx;
    const elt_vector_type& data    = m_x;

    // Compute the number of non-zero columns in each row
    idx_vector_type work (n_rows());
    std::for_each (row_idx.begin(), row_idx.end(), [&work](idx_type row) {work[row] ++;});
    
  
    idx_vector_type& b_col_idx = result.m_col_idx;
    idx_vector_type& b_row_idx = result.m_row_idx;
    elt_vector_type& b_data    = result.m_x;
    

    // Compute the column offsets for the transpose
    b_col_idx.resize (result.n_cols () + 1);
    b_col_idx[0] = 0;
    std::partial_sum (work.begin(), work.end(), b_col_idx.begin() + 1);
    std::copy (b_col_idx.begin(), b_col_idx.end() - 1, work.begin());
    
    // Compute the row indices and data
    b_row_idx.resize(nnz());
    b_data.resize   (nnz());
    
    // Iterate over the columns of the original matrix and update the
    // row indices of the transpose
    
    for (idx_type col = 0; col < n_cols(); col ++)
    {
        for (idx_type offset = col_idx[col]; offset < col_idx[col+1]; offset ++)
        {
            idx_type  row    = row_idx[offset];
            idx_type& index  = work[row];
            b_row_idx[index] = col;
            b_data   [index] = data[offset];
            index ++;
        }        
    }
}

template<class idx_type, class el_type>
csc_matrix<idx_type, el_type> csc_matrix<idx_type, el_type> :: transpose() const
{
  csc_matrix<idx_type, el_type> result (n_cols(), n_rows(), nnz());
  transpose(result);
  return result;
}

#endif // _CSC_MATRIX_TRANSPOSE_H_
