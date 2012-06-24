// -*- mode: c++ -*-
#ifndef _CSC_MATRIX_DECLARATIONS_H_
#define _CSC_MATRIX_DECLARATIONS_H_

#include "abstract_sparse_matrix.h"

/*! \brief A matrix using the compressed sparse column format.

	In compressed sparse column storage, the col_idx array is of size N + 1. col_idx[j] gives the starting position of the first non-zero element in column j. Hence col_idx[j+1] - col_idx[j] gives the total number of non-zero values in column j and therefore, col_idx[n_cols] gives the total number of non-zero elements in the matrix. 
	
	row_idx[j] and m_x[j] are arrays of size n_nzs, so col_idx[n_cols] == row_idx.size()
*/
	
template<class idx_type, class el_type> 
class csc_matrix : public abstract_sparse_matrix<idx_type, el_type>
{
public:
    typedef csc_matrix<idx_type, el_type> csc_t;

    using abstract_sparse_matrix<idx_type, el_type>::m_row_idx;
    using abstract_sparse_matrix<idx_type, el_type>::m_col_idx;
    using abstract_sparse_matrix<idx_type, el_type>::m_x;
    using abstract_sparse_matrix<idx_type, el_type>::m_n_rows;
    using abstract_sparse_matrix<idx_type, el_type>::m_n_cols;
    using abstract_sparse_matrix<idx_type, el_type>::n_rows;
    using abstract_sparse_matrix<idx_type, el_type>::n_cols;
    using abstract_sparse_matrix<idx_type, el_type>::shape;


    typedef typename abstract_sparse_matrix<idx_type, el_type>::idx_vector_type idx_vector_type;
    typedef typename abstract_sparse_matrix<idx_type, el_type>::elt_vector_type elt_vector_type;

    
private:
    static const idx_type npos;

public:
	
    csc_matrix (idx_type n_rows = 0, idx_type n_cols = 0, idx_type nz_max = 0): 
        abstract_sparse_matrix<idx_type, el_type> (n_rows, n_cols, nz_max) 
    {
		m_row_idx.reserve(nz_max);
        m_col_idx.resize (n_cols + 1);
    }
	
	/*! Returns number of nonzeros in the matrix. */
    virtual idx_type nnz() const
    {
        return m_row_idx.size();
    }  
	
	//----Matrix referencing/filling----//
	
	/*! Finds the (i,j)th coefficient of the matrix. Uses binary search to find an element.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\return The (i,j)th element of the matrix. 
	*/
	virtual el_type coeff(const idx_type& i, const idx_type& j) const {
		typename idx_vector_type::const_iterator low = lower_bound(m_row_idx.begin() + m_col_idx[j], m_row_idx.begin() + m_col_idx[j+1], i);
		return (*low == i) ? *(m_x.begin() + std::distance(m_row_idx.begin(), low)) : 0;
	}
	
	/*! Resizes the matrix. For use in preallocating space before factorization begins.
		\param n_rows the number of rows in the resized matrix.
		\param n_cols the number of cols in the resized matrix.
		\param n_nnzs the number of non-zeros expected in the matrix. 
	*/
	void resize(idx_type n_rows, idx_type n_cols, idx_type n_nzs)
	{
		m_n_rows = n_rows;
		m_n_cols = n_cols;
		m_col_idx.resize(n_cols + 1);
		m_row_idx.resize(n_nzs);
		m_x.resize(n_nzs);
	}
	
	//----Factorizations----//
	/*! Performs an LDL' factorization of this matrix. The factorization is performed in crout order and follows the algorithm outlined in "Crout versions of the ILU factorization with pivoting for sparse symmetric matrices" by Li and Saad (2005). Results are stored in L and D.
		\param L the L factor of this matrix.
		\param D the D factor of this matrix.
		\param lfil a parameter to control memory usage. Each column is guarannted to have fewer than lfil elements.
		\param tol a parameter to control agressiveness of dropping. In each column, elements less than tol*norm(column) are dropped.
	*/
	void ildl(csc_matrix<idx_type, el_type>& L, elt_vector_type& D, int lfil, double tol);
	
	//----IO Functions----//
	
    /*! \return A string reprepsentation of this matrix.
	*/
    std::string to_string () const;
	
	/*! \param filename the filename of the matrix to be loaded. Must be in matrix market format (.mtx).
	*/
	bool load(std::string filename);
	
	/*! \param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
	*/
	bool save(std::string filename);

};

#include "csc_matrix_ildl.h"
#include "csc_matrix_load.h"
#include "csc_matrix_save.h"
#include "csc_matrix_to_string.h"

template<class idx_type, class el_type>
const idx_type csc_matrix<idx_type, el_type>::npos = static_cast<idx_type>(-1);

#endif
