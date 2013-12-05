//-*- mode: c++ -*-
#ifndef _SKEW_BLOCK_DIAG_MATRIX_TO_STRING_H_
#define _SKEW_BLOCK_DIAG_MATRIX_TO_STRING_H_

#include <string>
#include <sstream>

template <class el_type>
std::string skew_block_diag_matrix<el_type> :: to_string() const
{
	std::ostringstream os;
	os << "Skew Block Diagonal Matrix (" << n_rows() << ", " << n_cols() << ", " << nnz() << ")" << std::endl;

	os << "Subdiagonal (col, val)  = " << "[";
	for (int i = 0; i < n_cols(); i=i+2)
		os << "(" << i << ", " << subdiag[i/2] << "), ";
	os << "]";

	return os.str();
}

#endif