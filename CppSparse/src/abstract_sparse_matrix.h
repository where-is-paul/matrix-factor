// -*- mode: c++ -*-
#ifndef _ABSTRACT_SPARSE_MATRIX_H_
#define _ABSTRACT_SPARSE_MATRIX_H_

#include <vector>
#include <string>

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
    virtual std::string to_string() const = 0;

    std::string __str__() const
    {
        return to_string();
    }

    std::string __repr__() const
    {
        return to_string();
    }

    virtual std::vector<el_type> full () const = 0;

    virtual ~abstract_sparse_matrix()
    {
    }
};

#endif // _ABSTRACT_SPARSE_MATRIX_H_
