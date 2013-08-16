// -*- mode: c++ -*-
#ifndef SYMMETRY_MATRIX_DECLARATIONS_H
#define SYMMETRY_MATRIX_DECLARATIONS_H


/*! \brief A list-of-lists (LIL) matrix in column oriented format.

	For convience, the matrix this class represents will be refered to as matrix A.
	In LIL-C format, each column of A (an n*n matrix) is stored as a separate vector. The nonzeros are stored in m_idx while the non-zeros are stored in m_x. Both m_x and m_idx are initialized to a list of n lists. m_idx and m_x are ordered dependent on each other, in that A(m_idx[k][j], k) = m_x[k][j].
	
*/
template <class el_type> 
class symmetry_matrix : public square_matrix<el_type>
{
public:
	
	//-------------- typedefs and inherited variables --------------//
	using square_matrix<el_type>::m_idx; // linux compiler gcc need these
	using square_matrix<el_type>::m_x;
	using square_matrix<el_type>::m_n_rows;
	using square_matrix<el_type>::m_n_cols;
	using square_matrix<el_type>::n_rows;
	using square_matrix<el_type>::n_cols;
	using square_matrix<el_type>::nnz;
	using square_matrix<el_type>::nnz_count;

	// using square_matrix<el_type>::idx_vector_type;
	// using square_matrix<el_type>::elt_vector_type;
	// using square_matrix<el_type>::idx_it;
	// using square_matrix<el_type>::elt_it;

	using square_matrix<el_type>::list;
	using square_matrix<el_type>::list_first;

	using square_matrix<el_type>::S;
	using square_matrix<el_type>::sign;

	typedef typename square_matrix<el_type>::idx_vector_type idx_vector_type; // linux compiler gcc need these
	typedef typename square_matrix<el_type>::elt_vector_type elt_vector_type;
	typedef typename square_matrix<el_type>::idx_it idx_it;
	typedef typename square_matrix<el_type>::elt_it elt_it;

	// typedef typename idx_vector_type::iterator idx_it;
	// typedef typename elt_vector_type::iterator elt_it;
	
public:
	
	/*! \brief Constructor for a column oriented list-of-lists (LIL) matrix. Space for both the values list and the indices list of the matrix is allocated here.
	*/
	symmetry_matrix(int n_rows = 0, int n_cols = 0) : square_matrix<el_type> (n_rows, n_cols)
	{
		sign = 1;
	}
	
	//----Matrix referencing/filling----//
	
	//----Factorizations----//
	/*! \brief Performs an LDL' factorization of this matrix. 
		
		The pivoted matrix P'AP will be stored in place of A. In addition, the L and D factors of P'AP will be stored in L and D (so that P'AP = LDL'). The factorization is performed in crout order and follows the algorithm outlined in "Crout versions of the ILU factorization with pivoting for sparse symmetric matrices" by Li and Saad (2005).
	
		\param L the L factor of this matrix.
		\param D the D factor of this matrix.
		\param perm the current permutation of A.
		\param fill_factor a parameter to control memory usage. Each column is guaranteed to have fewer than fill_factor*(nnz(A)/n_col(A)) elements.
		\param tol a parameter to control aggressiveness of dropping. In each column, elements less than tol*norm(column) are dropped.
    \param pp_tol a parameter to control aggresiveness of pivoting. Allowable ranges are [0,inf). If the parameter is >= 1, Bunch-Kaufman pivoting will be done in full. If the parameter is 0, partial pivoting will be turned off and the first non-zero pivot under the diagonal will be used. Choices close to 0 increase locality in pivoting (pivots closer to the diagonal are used) while choices closer to 1 increase the stability of pivoting. Useful for situations where you care more about preserving the structure of the matrix rather than bounding the size of its elements.
	*/
	void ildl(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol, const double& pp_tol);
	
	/*! \brief Performs an LDL' factorization of this matrix. 
		
		The pivoted matrix P'AP will be stored in place of A. In addition, the L and D factors of P'AP will be stored in L and D (so that P'AP = LDL'). The factorization is performed in crout order and follows rook pivoting strategy.
	
		\param L the L factor of this matrix.
		\param D the D factor of this matrix.
		\param perm the current permutation of A.
		\param fill_factor a parameter to control memory usage. Each column is guaranteed to have fewer than fill_factor*(nnz(A)/n_col(A)) elements.
		\param tol a parameter to control agressiveness of dropping. In each column, elements less than tol*norm(column) are dropped.
	*/
	void ildlrp(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol);
	
	//------Helpers------//
	
	/*! \brief Updates A.list for iteration k.
		\param k current iteration index.	
	*/
	inline void advance_list(const int& k)
	{
		for (idx_it it = m_idx[k].begin(); it != m_idx[k].end(); it++)
		{
			if (*it == k) continue;
			this->ensure_invariant(*it, k, list[*it]); //make sure next element is good.
			list_first[*it]++; //invariant ensured.
		}
	}
	
	inline void update_single(const int& j, const el_type& l_ki, const el_type& d, std::vector<el_type>& work, std::vector<int>& curr_nnzs, ultriangular_matrix<el_type>& L, vector<bool>& in_set);
	inline void update(const int& r, std::vector<el_type>& work, std::vector<int>& curr_nnzs, ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, vector<bool>& in_set);
	inline void calculate(int& k, ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, elt_vector_type& col_i, elt_vector_type& col_r, idx_vector_type& col_i_nnzs, idx_vector_type& col_r_nnzs, bool& size_two_piv, const double& tol, int& lfil, vector<bool>& in_set, const double& stat_piv);

	/*!	\brief Finds the degree of a node.
		\param index the index of the node.
		\return the degree of the nodex indexed by index.
	*/
	inline int degree(const int& index) const
	{
		int deg = list[index].size() + m_idx[index].size();
		if (m_idx[index].size() > 0 && m_idx[index][0] == index) deg--;
		return deg;
	}
	
	//----IO helper Functions----//
	// overwrite the helper functions in the class square_matrix

	void put_header(std::string& header)
	{
		header = "%%MatrixMarket matrix coordinate real symmetric";
	}

	inline bool readline(char*& line, int& n_rows, int& n_cols, int& i, int& j, el_type& value)
	{
		char * end1, * end2;
		i = strtol(line, &end1, 10);
		j = strtol(end1, &end2, 10);
		
		double val;
		val = strtod(end2, NULL); //change right here if you want to read in complex or somethin
		value = val;
		
		i--; j--;
		return (i>=0 && j>=0 && i<n_rows&& j<n_cols);
	}
};

//------------------ include files for class functions -------------------//

#include "symmetry_matrix_update.h"
#include "symmetry_matrix_ildl.h"
#include "symmetry_matrix_ildlrp.h"
#include "symmetry_matrix_calculate.h"

#endif