//-*- mode: c++ -*-
#ifndef _BLOCK_DIAG_MATRIX_TO_STRING_H_
#define _BLOCK_DIAG_MATRIX_TO_STRING_H_

#include <string>
#include <sstream>

template <class el_type>
std::string block_diag_matrix<el_type> :: to_string() const
{
  std::ostringstream os;
  os << "Block Diagonal Matrix (" << n_rows() << ", " << n_cols() << ", " << nnz() << ")" << std::endl;
  
  os << "Main Diagonal Values          = " << main_diag << std::endl;
  os << "Off Diagonal (Column, Value)  = " << "[";
  for (int i = 0; i < n_cols(); i++) {
    if (block_size(i) == 2) {
		os << "(" << i << ", " << off_diag.find(i)->second << "), ";
		i++;
	}
  }
  os << "]";
  
  
  return os.str();
}

#endif
