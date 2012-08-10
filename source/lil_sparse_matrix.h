// -*- mode: c++ -*-
#ifndef _LIL_SPARSE_MATRIX_H_
#define _LIL_SPARSE_MATRIX_H_

#include <vector>
#include <string>
#include <fstream>
#include <limits>

#ifndef VECTOR_SHIFT
#define VECTOR_SHIFT
template<class Container>
std::ostream& operator<< (std::ostream& os, const Container& vec)
{
	os << "[";
	if (!vec.empty())
	{
		for (auto it = vec.begin(); it+1 != vec.end(); it++)
		{
			os << *it << ", ";
		}
		
		it++;
		os << *it;
	}
	os << "]";
	return os;
}
#endif

using std::vector;

/*! \brief The abstract parent of all sparse matrices */
template<class el_type>
class lil_sparse_matrix 
{

public:

	typedef vector<int> idx_vector_type;
	typedef vector<el_type>  elt_vector_type;
	
	/*! \brief Allows outputting the contents of the matrix via << operators. */
	friend std::ostream & operator<<(std::ostream& os, const lil_sparse_matrix& A) 
	{
		os << A.to_string();
		return os;
	};

	/*! Number of rows/cols in the matrix */
	int m_n_rows, m_n_cols, nnz_count;
	el_type eps;

	/*! The row/col indices. The way m_idx is used depends on whether the matrix is in LIL-C or LIL-R. */
	vector<idx_vector_type> m_idx;
	/*! The values of the nonzeros in the matrix. */
	vector<elt_vector_type> m_x;

	/*! \brief Default constructor for an abstract matrix. This constructor will be extended by base classes depending on the representation of the matrix (LIL-C or LIL-R). */
	lil_sparse_matrix (int n_rows, int n_cols) : m_n_rows(n_rows), m_n_cols (n_cols)
	{
		nnz_count = 0;
		eps = std::numeric_limits<el_type>::epsilon();
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
	
	/*! \brief Returns A_ij (zero-indexed). This function should be extended by subclasses as it is dependent on the matrix storage type.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\return The (i,j)th element of the matrix. */
	virtual el_type coeff(const int& i, const int& j) const = 0;
	
	/*! \return A string reprepsentation of this matrix.
	*/
	virtual std::string to_string() const = 0;

	/*! Generic class destructor. */
	virtual ~lil_sparse_matrix()
	{
	}
};

#endif // _LIL_SPARSE_MATRIX_H_
