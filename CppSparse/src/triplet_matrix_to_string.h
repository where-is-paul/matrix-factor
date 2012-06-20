//-*-mode:c++-*-
#ifndef _TRIPLET_MATRIX_TO_STRING_H_
#define _TRIPLET_MATRIX_TO_STRING_H_

#include <string>
#include <sstream>

template <class idx_type, class el_type>
std::string triplet_matrix<idx_type, el_type> :: to_string () const
{
  std::ostringstream os;
  os << "Triplet Sparse Matrix (" << m_n_rows << ", " << m_n_cols << ", " << this->nz_max() << ")" << std::endl;

  for (size_t i = 0; i < m_row_idx.size(); i++) 
  {
    os << "\t (" << m_row_idx[i] << ", " << m_col_idx[i] << ") = " << m_x[i] << std::endl;
  }
  
  return os.str();
}


#endif
