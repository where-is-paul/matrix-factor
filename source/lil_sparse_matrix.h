// -*- mode: c++ -*-
#ifndef LIL_SPARSE_MATRIX_H
#define LIL_SPARSE_MATRIX_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>

#include "lil_sparse_matrix_helpers.h"

/*! \brief The abstract parent of all sparse matrices */
template<class el_type>
class lil_sparse_matrix 
{

public:

	typedef vector<int> idx_vector_type;
	typedef vector<el_type> elt_vector_type;
	typedef idx_vector_type::iterator idx_it;
	typedef typename elt_vector_type::iterator elt_it;

	/*! \brief Allows outputting the contents of the matrix via << operators. */
	friend std::ostream & operator<<(std::ostream& os, const lil_sparse_matrix& A) 
	{
		os << A.to_string();
		return os;
	};

	int m_n_rows;///<Number of rows in the matrix.
	int	m_n_cols;///<Number of cols in the matrix.
	int nnz_count;///<Number of nonzeros in the matrix.

	vector<idx_vector_type> m_idx;///<The row/col indices. The way m_idx is used depends on whether the matrix is in LIL-C or LIL-R.
	vector<elt_vector_type> m_x;///<The values of the nonzeros in the matrix.
	vector<idx_vector_type> list;	///<A list of linked lists that gives the non-zero elements in each row of A. Since at any time we may swap between two rows, we require linked lists for each row of A.
	
	/*! \brief Default constructor for an abstract matrix. This constructor will be extended by base classes depending on the representation of the matrix (LIL-C or LIL-R). */
	lil_sparse_matrix (int n_rows, int n_cols) : m_n_rows(n_rows), m_n_cols (n_cols)
	{
		nnz_count = 0;
	}
	
	/*! \return Number of rows in the matrix. */
	int n_rows() const
	{
		return m_n_rows;
	}

	/*! \return Number of cols in the matrix. */
	int n_cols() const
	{
		return m_n_cols;
	}

	/*! \return Number of nonzeros in the matrix. */
	int nnz() const
	{
		return nnz_count;
	};

	/*! \brief Finds the (i,j)th coefficient of the matrix.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\param offset an optional search offset for use in linear search (start at offset instead of 0).
		\return The (i,j)th element of the matrix.
	*/
	inline el_type coeff(const int& i, const int& j, int offset = 0) const
	{	
		for (unsigned int k = offset, end = m_idx[j].size(); k < end; k++)
			if (m_idx[j][k] == i)
				return m_x[j][k];
		
		return 0;
	}

	/*! \brief Resizes the matrix. For use in preallocating space before factorization begins.
		\param n_rows the number of rows in the resized matrix.
		\param n_cols the number of cols in the resized matrix.
	*/
	void resize(int n_rows, int n_cols)
	{
		m_n_rows = n_rows;
		m_n_cols = n_cols;
		m_idx.resize(n_cols);
		m_x.resize(n_cols);
		list.resize(n_cols);
	}

	/*! \return A string reprepsentation of this matrix.
	*/
	virtual std::string to_string() const = 0;

	/*! Generic class destructor. */
	virtual ~lil_sparse_matrix()
	{
	}
};

#endif