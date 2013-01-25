//-*- mode: c++ -*-
#ifndef _BLOCK_DIAG_MATRIX_TO_STRING_H_
#define _BLOCK_DIAG_MATRIX_TO_STRING_H_

#include <string>
#include <sstream>

template <class el_type>
std::string block_diag_matrix<el_type> :: to_string() const
{
	bool first = true;
	std::ostringstream os;
	os << "Block Diagonal Matrix (" << n_rows() << ", " << n_cols() << ", " << nnz() << ")" << std::endl;
  
	os << "Main Diagonal Values     = [" << main_diag << "]" << std::endl;
	os << "Off Diagonal (col, val)  = [";
	for (int i = 0; i < n_cols(); i++)
	{
		if (block_size(i) == 2)
		{
			if (!first)
				os << ", (" << i << ", " << off_diag.find(i)->second << ")";
			else
			{
				first = false;
				os << "(" << i << ", " << off_diag.find(i)->second << ")";
			}
			i++;
		}
	}
	os << "]" << std::endl;

	return os.str();
}

#endif