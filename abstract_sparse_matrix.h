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

template<class idx_type, class el_type>
class abstract_sparse_matrix 
{
public:
    typedef idx_type index_type;
    typedef el_type  element_type;
  
    typedef std::vector<idx_type> idx_vector_type;
    typedef std::vector<el_type>  elt_vector_type;
	
	friend std::ostream & operator<<(std::ostream& os, const abstract_sparse_matrix& A) 
	{
		os << A.to_string();
		return os;
	};

protected:
  
    idx_type m_n_rows, m_n_cols;

    idx_vector_type m_col_idx, m_row_idx;
    elt_vector_type m_x;

    abstract_sparse_matrix (idx_type n_rows, idx_type n_cols, idx_type nz_max = 0) : m_n_rows(n_rows), m_n_cols (n_cols)
    {
        m_row_idx.reserve(nz_max);
        m_x.reserve(nz_max);
    }

public:
    /// 
    /// 
    ///
    const idx_type n_rows() const
    {
        return m_n_rows;
    }

    const idx_type n_cols() const
    {
        return m_n_cols;
    }

    const std::pair<idx_type, idx_type> shape() const
    {
        return std::make_pair (m_n_rows, m_n_cols);
    }
    
    const idx_type nz_max() const
    {
        return m_x.capacity();
    }

    const idx_vector_type& row_idx () const
    {
        return m_row_idx;
    }

    const idx_vector_type& col_idx () const
    {
        return m_col_idx;
    }

    const elt_vector_type& nz_vals () const
    {
        return m_x;
    }

    virtual idx_type nnz() const = 0;
    virtual el_type coeff(const idx_type& i, const idx_type& j) const = 0;
    virtual std::string to_string() const = 0;


    virtual ~abstract_sparse_matrix()
    {
    }
};

#endif // _ABSTRACT_SPARSE_MATRIX_H_
