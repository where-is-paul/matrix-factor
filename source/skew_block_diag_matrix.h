// -*- mode: c++ -*-
#ifndef _SKEW_BLOCK_DIAG_MATRIX_H_
#define _SKEW_BLOCK_DIAG_MATRIX_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>

#ifndef DEBUG
#define DEBUG
template<class el_type>
std::ostream& operator<< (std::ostream& os, const std::vector<el_type>& vec)
{
	os << "[";
	if (!vec.empty())
	{
		for (typename std::vector<el_type>::size_type index = 0; index < vec.size() - 1; index ++)
		{
			os << vec[index] << ", ";
		}

		os << vec[vec.size()-1];
	}
	os << "]";
	return os;
}
#endif


/*! \brief A quick implementation of a diagonal matrix with 1x1 and 2x2 blocks. 
*/
template<class el_type>
class skew_block_diag_matrix
{
public:

	typedef std::vector<el_type>  elt_vector_type;
	
	/*! Allows outputting the contents of the matrix via << operators. */
	friend std::ostream& operator<< (std::ostream& os, const skew_block_diag_matrix& D) 
	{
		os << D.to_string();
		return os;
	};
	
	int m_n_size;///<Dimension of the matrix.
	int nnz_count;///<Number of non-zeros in the matrix.
	elt_vector_type subdiag;///<Stores subdiagonal elements.
	
	/*!	\brief Constructor for diagonal class. Initializes a 0x0 matrix when given no arguments.
	*/
	skew_block_diag_matrix (int n_rows = 0, int n_cols = 0) : m_n_size(n_rows)
	{
		assert(n_rows == n_cols);
		nnz_count = n_rows / 2;
		subdiag.resize(n_rows/2);
	}
	
	/*!	\brief Resizes this matrix to an n*n matrix.
	*/
	void resize(int n)
	{
		m_n_size = n;
		subdiag.resize(n/2);
		nnz_count = n / 2;
	}
	
	/*! \return Number of rows in the matrix. */
	int n_rows() const
	{
		return m_n_size;
	}

	/*! \return Number of cols in the matrix. */
	int n_cols() const
	{
		return m_n_size;
	}

	/*! \return Number of nonzeros in the matrix. */
	int nnz() const 
	{
		return nnz_count;
	};
	
	/*!	\param i the index of the element.
		\return The D(i+1,i)th element if i is even, or D(i, i-1) if i is odd
	*/
	el_type& operator[](int i)
	{
		return subdiag.at(i / 2);
	}
	
	/*! \return A string reprepsentation of this matrix.
	*/
	std::string to_string() const;
	
	/*! \param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
		\return True if the save succeeded, false otherwise.
	*/
	bool save(std::string filename) const;
	
	/*! \brief Generic class destructor.
	*/
	~skew_block_diag_matrix()
	{
	}
};

#include "skew_block_diag_matrix_to_string.h"
#include "skew_block_diag_matrix_save.h"

#endif 