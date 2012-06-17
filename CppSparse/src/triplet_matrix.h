//-*- mode: c++ -*-
#ifndef _TRIPLET_MATRIX_H_
#define _TRIPLET_MATRIX_H_
#include "abstract_sparse_matrix.h"

template<class idx_type, class el_type>
class triplet_matrix : public abstract_sparse_matrix<idx_type, el_type>
{
    using abstract_sparse_matrix<idx_type, el_type>::m_row_idx;
    using abstract_sparse_matrix<idx_type, el_type>::m_col_idx;
    using abstract_sparse_matrix<idx_type, el_type>::m_x;
    using abstract_sparse_matrix<idx_type, el_type>::m_n_rows;
    using abstract_sparse_matrix<idx_type, el_type>::m_n_cols;

public:
    using abstract_sparse_matrix<idx_type, el_type>::n_rows;
    using abstract_sparse_matrix<idx_type, el_type>::n_cols;

    triplet_matrix (idx_type rows = 0, idx_type cols = 0, idx_type nz_max = 0) : abstract_sparse_matrix<idx_type, el_type> (rows, cols, nz_max)
    {
        m_col_idx.reserve (nz_max);
    }
    
    explicit triplet_matrix (std::istream&);
    explicit triplet_matrix (const std::string& );
    
    void push_back (idx_type row, idx_type col, el_type value)
    {
        m_row_idx.push_back(row);
        m_col_idx.push_back(col);
        m_x.push_back(value);
        
        m_n_rows = std::max(m_n_rows, row + 1);
        m_n_cols = std::max(m_n_cols, col + 1);
    }

    void push_back (const std::vector<idx_type>& rows, const std::vector<idx_type>& cols, const std::vector<el_type>& values)
    {
        idx_type k = static_cast<idx_type> (rows.size());
        std::cout << "k = " << k << std::endl;
        if (!(k == cols.size() && (k * k) == values.size()))
        {
            throw std::logic_error(__FUNCTION__);
        }

        for (idx_type j = 0; j < k; j ++)
        {
            for (idx_type i = 0; i < k; i ++)
            {
                std::cout << "(" << i << ", " << j << ")" << std::endl;
                push_back(rows[i], cols[j], values[i + k * j]);
            }
        }
    }
    
    virtual idx_type nnz () const 
    {
        return m_x.size();
    }

    void gaxpy (const el_type* x, el_type* y) const;

    void gaxpy (const std::vector<el_type>& x, std::vector<el_type>& y) const;

    virtual std::string to_string() const;    
    
    virtual std::vector<el_type> full () const;
private:
    void load(std::istream&);
};

#include "triplet_matrix_load.h"
#include "triplet_matrix_to_string.h"
#include "triplet_matrix_full.h"
#include "triplet_matrix_gaxpy.h"

#endif
