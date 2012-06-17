//-*- mode: c++ -*-
#ifndef _CSC_MATRIX_SUM_DUPLICATES_H_
#define _CSC_MATRIX_SUM_DUPLICATES_H_

template<class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: sum_duplicates ()
{
  // Initially the row offsets are unassigned
  idx_vector_type new_offsets (n_rows(), npos);

  idx_type new_offset = 0;

  for (idx_type col = 0; col < n_cols(); col ++)
  {
      idx_type col_start = new_offset;

      for (idx_type old_offset = m_col_idx [col]; old_offset < m_col_idx[col+1]; old_offset ++)
      {
          idx_type row  = m_row_idx[old_offset];
          el_type  data = m_x      [old_offset];
          
          if (new_offsets[row] == npos || new_offsets[row] < col_start)
          {
              m_row_idx  [new_offset] = row;
              m_x        [new_offset] = data;
              new_offsets[row]        = new_offset;
              new_offset ++;
          }
          else
          {
              m_x[new_offsets[row]] += data;
          }
      }
      m_col_idx[col] = col_start;
  }
  m_col_idx.back() = new_offset;
  m_row_idx.resize(new_offset);
  m_x.resize(new_offset);
}

#endif //_CSC_MATRIX_SUM_DUPLICATES_H_
