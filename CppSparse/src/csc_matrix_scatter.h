//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_SCATTER_H_
#define _CSC_MATRIX_SCATTER_H_

/// scatter computes:
/// x += A(:, k) * beta
/// This is in turn used to compute each column of C as
/// C(:, j) += A(:, k) * b(k, j) 

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: scatter (csc_matrix<idx_type, el_type>& result, idx_type column, el_type beta, 
					       idx_vector_type& work, elt_vector_type& x, idx_type mark, idx_type& nnz_result) const
{
    for (idx_type offset = m_col_idx[column]; offset < m_col_idx[column+1]; offset ++)
    {
	idx_type row = m_row_idx[offset];
	el_type  val = m_x[offset];
	if (work[row] < mark)
	{
	    work[row] = mark;
	    result.m_row_idx.push_back(row);
            result.m_x.push_back(0);
	    nnz_result ++;
	}
	else
	{
	    x[row] += val * beta;
	}
    }
}

#endif // _CSC_MATRIX_SCATTER_H_
