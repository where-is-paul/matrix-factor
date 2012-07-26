// -*- mode: c++ -*-
#ifndef _LILC_MATRIX_DECLARATIONS_H_
#define _LILC_MATRIX_DECLARATIONS_H_

#include "lil_sparse_matrix.h"
#include <deque>
#include <algorithm>
#include <cmath>

/*! \brief A matrix using the compressed sparse column format.

	In compressed sparse column storage, the col_idx array is of size N + 1. col_idx[j] gives the starting position of the first non-zero element in column j. Hence col_idx[j+1] - col_idx[j] gives the total number of non-zero values in column j and therefore, col_idx[n_cols] gives the total number of non-zero elements in the matrix. 
	
	row_idx[j] and m_x[j] are arrays of size n_nzs, so col_idx[n_cols] == row_idx.size()
*/

template<class el_type> 
class lilc_matrix : public lil_sparse_matrix<el_type>
{
public:

	using lil_sparse_matrix<el_type>::m_idx;
	using lil_sparse_matrix<el_type>::m_x;
	using lil_sparse_matrix<el_type>::m_n_rows;
	using lil_sparse_matrix<el_type>::m_n_cols;
	using lil_sparse_matrix<el_type>::n_rows;
	using lil_sparse_matrix<el_type>::n_cols;
	using lil_sparse_matrix<el_type>::nnz;
	using lil_sparse_matrix<el_type>::nnz_count;


	typedef typename lil_sparse_matrix<el_type>::idx_vector_type idx_vector_type;
	typedef typename lil_sparse_matrix<el_type>::elt_vector_type elt_vector_type;
	
	typedef typename idx_vector_type::iterator idx_it;
	typedef typename elt_vector_type::iterator elt_it;
	typedef std::deque< int > ::iterator list_it;

	std::vector< std::deque< int > > list;
	std::vector<int> first;
	
	block_diag_matrix<el_type> S;
	
public:
	
	lilc_matrix (int n_rows = 0, int n_cols = 0): 
	lil_sparse_matrix<el_type> (n_rows, n_cols) 
	{
		m_x.reserve(n_cols);
		m_idx.reserve(n_cols);
	}
	
	//----Matrix referencing/filling----//
	
	/*! Finds the (i,j)th coefficient of the matrix.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\return The (i,j)th element of the matrix. 
	*/
	virtual el_type coeff(const int& i, const int& j) const 
	{	
		//invariant: first elem in each col of a is the diagonal elem if it exists.
		if (i == j) {
			if (m_idx[j].size() == 0) return 0;
			return (m_idx[j][0] == i ? m_x[j][0] : 0);
		}
		
		for (unsigned int k = 0; k < m_idx[j].size(); k++) {
			if (m_idx[j][k] == i) return m_x[j][k];
		}
		
		return 0;
	}
	
	/*! Finds the (i,j)th coefficient of the matrix.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\return True if (i,j)th element is nonzero, false otherwise. 
	*/
	bool coeffRef(const int& i, const int& j, std::pair<idx_it, elt_it>& its)
	{	
		for (unsigned int k = 0; k < m_idx[j].size(); k++) {
			if (m_idx[j][k] == i) {
				its = make_pair(m_idx[j].begin() + k, m_x[j].begin() + k);
				return true;
			}
		}
		
		its = make_pair(m_idx[j].end(), m_x[j].end());
		return false;
	}
	
	/*! Resizes the matrix. For use in preallocating space before factorization begins.
		\param n_rows the number of rows in the resized matrix.
		\param n_cols the number of cols in the resized matrix.
		\param n_nzs the number of non-zeros expected in the matrix. 
	*/
	void resize(int n_rows, int n_cols)
	{
		m_n_rows = n_rows;
		m_n_cols = n_cols;
		
		m_x.resize(n_cols);
		m_idx.resize(n_cols);
		
		first.resize(n_cols, 1);
		list.resize(n_cols);
		
		S.resize(n_cols);
	}
	//-----Reorderings/Rescalings------//
	void sym_perm(vector<int>& perm);
	
	/*!the symmetric  matrix A of order n is equilibrated  and  the symmetric  equilibrated  matrix  SAS  is stored  in the  lower triangular  part  of A,  where  S-1 = diag (S[1],..., S[n]);
	*/
	void sym_equil();
	
	//----Factorizations----//
	/*! Performs an LDL' factorization of this matrix. The factorization is performed in crout order and follows the algorithm outlined in "Crout versions of the ILU factorization with pivoting for sparse symmetric matrices" by Li and Saad (2005). Results are stored in L and D.
		\param L the L factor of this matrix.
		\param D the D factor of this matrix.
		\param lfil a parameter to control memory usage. Each column is guarannted to have fewer than lfil elements.
		\param tol a parameter to control agressiveness of dropping. In each column, elements less than tol*norm(column) are dropped.
	*/
	void ildl(lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, double fill_factor, double tol);
	
	//------Helpers------//
	inline void pivot(vector<idx_it>& swapk, vector<idx_it>& swapr, vector<list_it>& swapk_, vector<list_it>& swapr_, idx_vector_type& all_swaps, vector<bool>& in_set, elt_vector_type& col_k, idx_vector_type& col_k_nnzs, elt_vector_type& col_r, idx_vector_type& col_r_nnzs, lilc_matrix<el_type>& L, const int& k, const int& r);
	
	template <class Container>
	inline void ensure_invariant(const int& j, const int& k, Container& con, int offset, bool update_list = false) {
		if ((offset >= (int) con.size()) || con.empty() || con[offset] == k) return;
		
		int i, min(offset);
		for (i = offset; i < (int) con.size(); i++) {
			if (con[i] == k) {
				min = i; 
				break;
			} else if ( con[i] < con[min] ) {
				min = i;
			}
		}
		
		if (update_list)
			std::swap(con[offset], con[min]);
		else {
			std::swap(con[first[j]], con[min]);
			std::swap(m_x[j][first[j]], m_x[j][min]);
		}
	}
	
	inline void advance_first(const int& k) {
		int offset;
		for (auto it = list[k].begin(); it != list[k].end(); it++) {	
			//ensure invariant (perhaps not needed here since the invariant is always ensured during pivoting).
			offset = first[*it];

			if (offset >= (int) m_idx[*it].size()) continue;
			
			ensure_invariant(*it, k, m_idx[*it], offset);
			
			if (m_idx[*it][offset] == k)
				first[*it]++;
		}
	}
	
	inline void advance_list(const int& k) {
			int offset;
			for (auto it = m_idx[k].begin(); it != m_idx[k].end(); it++) {
				offset = first[*it];
				

				if (offset >= (int) list[*it].size()) continue;
				
				ensure_invariant(*it, k, list[*it], offset, true);
				
				if (list[*it][offset] == k)
					first[*it]++;
			}
					
	}
	
	//----IO Functions----//
	
	/*! \return A string reprepsentation of this matrix.
	*/
	std::string to_string () const;
	
	/*! \param filename the filename of the matrix to be loaded. Must be in matrix market format (.mtx).
	*/
	bool load(std::string filename);
	
	/*! \param filename the filename of the matrix to be saved. All matrices saved are in matrix market format (.mtx).
	*/
	bool save(std::string filename, bool sym);

};

#include "lilc_matrix_sym_perm.h"
#include "lilc_matrix_sym_equil.h"
#include "lilc_matrix_ildl.h"
#include "lilc_matrix_ildl_helpers.h"
#include "lilc_matrix_pivot.h"
#include "lilc_matrix_load.h"
#include "lilc_matrix_save.h"
#include "lilc_matrix_to_string.h"

#endif
