// -*- mode: c++ -*-
#ifndef _ABSTRACT_SPARSE_MATRIX_H_
#define _ABSTRACT_SPARSE_MATRIX_H_

#include <vector>
#include <string>
#include <fstream>

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

/*! \brief The abstract parent of all sparse matrices */
template<class idx_type, class el_type>
class abstract_sparse_matrix 
{
public:
	/*! The index type of the matrix. */
    typedef idx_type index_type;
	/*! The element type of the matrix (double, complex, etc). */
    typedef el_type  element_type;
  
    typedef std::vector<idx_type> idx_vector_type;
    typedef std::vector<el_type>  elt_vector_type;
	
	/*! Allows outputting the contents of the matrix via << operators. */
	friend std::ostream & operator<<(std::ostream& os, const abstract_sparse_matrix& A) 
	{
		os << A.to_string();
		return os;
	};

protected:
  
	/*! Number of rows/cols in the matrix */
    idx_type m_n_rows, m_n_cols;

	/*! The row/col indices. The way m_col_idx and m_row_idx are used depends on whether the matrix is in CSC, CSR, or Triplet form. */
    idx_vector_type m_col_idx, m_row_idx;
	/*! The values of the nonzeros in the matrix. */
    elt_vector_type m_x;

	/*! Default constructor for an abstract matrix. This constructor will be extended by base classes depending on the representation of the matrix (CSC, CSR, or Triplet). */
    abstract_sparse_matrix (idx_type n_rows, idx_type n_cols, idx_type nz_max = 0) : m_n_rows(n_rows), m_n_cols (n_cols)
    {
        m_x.reserve(nz_max);
    }

public:
    
	/*! \return Number of rows in the matrix. */
    const idx_type n_rows() const
    {
        return m_n_rows;
    }

	/*! \return Number of cols in the matrix. */
    const idx_type n_cols() const
    {
        return m_n_cols;
    }

	/*! \return Number of rows/cols in the matrix. */
    const std::pair<idx_type, idx_type> shape() const
    {
        return std::make_pair (m_n_rows, m_n_cols);
    }
    
	/*! \return The maximum number of nonzeros this matrix may hold without resizing. */
    const idx_type nz_max() const
    {
        return m_x.capacity();
    }

	/*! \return All row indices in the matrix. */
    const idx_vector_type& row_idx () const
    {
        return m_row_idx;
    }

	/*! \return All col indices in the matrix. */
    const idx_vector_type& col_idx () const
    {
        return m_col_idx;
    }

	/*! \return All nonzeros in the matrix. */
    const elt_vector_type& nz_vals () const
    {
        return m_x;
    }

	/*! \return Number of nonzeros in the matrix. */
    virtual idx_type nnz() const = 0;
	
	/*! Returns A_ij (zero-indexed). This function should be extended by subclasses as it is dependent on the matrix storage type.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\return The (i,j)th element of the matrix. */
    virtual el_type coeff(const idx_type& i, const idx_type& j) const = 0;
	
	/*! Returns a string representation of the matrix. */
    virtual std::string to_string() const = 0;

	/*! Generic class destructor. */
    virtual ~abstract_sparse_matrix()
    {
    }
};

#endif // _ABSTRACT_SPARSE_MATRIX_H_
