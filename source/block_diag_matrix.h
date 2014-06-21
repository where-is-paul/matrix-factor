// -*- mode: c++ -*-
#ifndef _BLOCK_DIAG_MATRIX_H_
#define _BLOCK_DIAG_MATRIX_H_

#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>

#ifdef SYM_ILDL_DEBUG
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
class block_diag_matrix
{
public:

	typedef std::unordered_map<int, el_type> int_elt_map;
	typedef std::vector<el_type>  elt_vector_type;
	
	/*! Allows outputting the contents of the matrix via << operators. */
	friend std::ostream & operator<<(std::ostream& os, const block_diag_matrix& D) 
	{
		os << D.to_string();
		return os;
	};
	
	int m_n_size;///<Dimension of the matrix.
	int nnz_count;///<Number of non-zeros in the matrix.
	elt_vector_type main_diag;///<Stores main diagonal elements.
	int_elt_map off_diag;///<Stores off-diagonal elements of 2x2 pivots.
	
	/*!	\brief Constructor for diagonal class. Initializes a 0x0 matrix when given no arguments.
	*/
	block_diag_matrix (int n_rows = 0, int n_cols = 0) : m_n_size(n_rows)
	{
		assert(n_rows == n_cols);
		nnz_count = n_rows;
		main_diag.resize(n_rows);
	}
	
	/*!	\brief Resizes this matrix to an n*n matrix with default_value on the main diagonal.
	*/
	void resize(int n, el_type default_value) {
		m_n_size = n;
		main_diag.resize(n, default_value);
		nnz_count = n;
	}
	
	/*!	\brief Resizes this matrix to an n*n matrix.
	*/
	void resize(int n) {
		resize(n, 0);
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
		\return The D(i,i)th element.
	*/
	el_type& operator[](int i) {
		return main_diag.at(i);
	}
	
	/*!	\param i the index of the element.
		\return The D(i+1,i)th element.
	*/
	el_type& off_diagonal(int i) {
		if (!off_diag.count(i)) {
			off_diag.insert(std::make_pair(i, 0));
			nnz_count++;
		}
		
		return off_diag[i];
	}
	
	/*!	\param i the index of the element.
		\return 2 if there is a diagonal pivot at D(i,i) and D(i+1,i+1).
				-2 if there is a diagonal pivot at D(i-1,i-1) and D(i,i).
				1 if the pivot is only a 1x1 block.
	*/
	int block_size(int i) const {
		if (off_diag.count(i)) {
			return 2;
		} else if (off_diag.count(i-1)) {
			return -2;
		} else {
			return 1;
		}
	}
	
	/*! \return A string reprepsentation of this matrix.
	*/
	std::string to_string () const;
	
	/*! \param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
		\return True if the save succeeded, false otherwise.
	*/
	bool save(std::string filename) const;
	
	/*! \brief Generic class destructor.
	*/
	~block_diag_matrix()
	{
	}
};

#include "block_diag_matrix_to_string.h"
#include "block_diag_matrix_save.h"

#endif 