// -*- mode: c++ -*-
#ifndef _LIL_SPARSE_MATRIX_H_
#define _LIL_SPARSE_MATRIX_H_

#include <vector>
#include <string>
#include <fstream>
#include <limits>

template<class el_type>
bool save(const std::vector<el_type>& vec, std::string filename) {
	std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary);
	if(!out)
	return false;

	out.flags(std::ios_base::scientific);
	out.precision(12);
	std::string header = "%%MatrixMarket matrix coordinate real general";; 

	out << header << std::endl; 
	out << vec.size() << " " << vec.size() << " " << vec.size() << "\n";

	for(int i = 0; i < (int) vec.size(); i++)
	out << i+1 << " " << i+1 << " " << vec[i] << "\n";
	
	out.close();
	return true;
}

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
	/*! The element type of the matrix (double, complex, etc). */
	typedef el_type  element_type;

	typedef vector<int> idx_vector_type;
	typedef vector<el_type>  elt_vector_type;
	
	/*! Allows outputting the contents of the matrix via << operators. */
	friend std::ostream & operator<<(std::ostream& os, const lil_sparse_matrix& A) 
	{
		os << A.to_string();
		return os;
	};

protected:
	
	/*! Number of rows/cols in the matrix */
	int m_n_rows, m_n_cols, nnz_count;
	el_type eps;

	/*! The row/col indices. The way m_col_idx and m_row_idx are used depends on whether the matrix is in CSC, CSR, or Triplet form. */
	vector<idx_vector_type> m_idx;
	/*! The values of the nonzeros in the matrix. */
	vector<elt_vector_type> m_x;

	/*! Default constructor for an abstract matrix. This constructor will be extended by base classes depending on the representation of the matrix (CSC, CSR, or Triplet). */
	lil_sparse_matrix (int n_rows, int n_cols) : m_n_rows(n_rows), m_n_cols (n_cols)
	{
		nnz_count = 0;
		eps = std::numeric_limits<el_type>::epsilon();
	}

public:
	
	/*! \return Number of rows in the matrix. */
	const int n_rows() const
	{
		return m_n_rows;
	}

	/*! \return Number of cols in the matrix. */
	const int n_cols() const
	{
		return m_n_cols;
	}

	/*! \return Number of nonzeros in the matrix. */
	virtual int nnz() const 
	{
		return nnz_count;
	};
	
	/*! Returns A_ij (zero-indexed). This function should be extended by subclasses as it is dependent on the matrix storage type.
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
