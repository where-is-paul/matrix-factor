// -*- mode: c++ -*-
#ifndef SKEW_SYMMETRY_MATRIX_DECLARATIONS_H
#define SKEW_SYMMETRY_MATRIX_DECLARATIONS_H

/*! \brief A list-of-lists (LIL) matrix in column oriented format.

	For convience, the matrix this class represents will be refered to as matrix A.
	In LIL-C format, each column of A (an n*n matrix) is stored as a separate vector. The nonzeros are stored in m_idx while the non-zeros are stored in m_x. Both m_x and m_idx are initialized to a list of n lists. m_idx and m_x are ordered dependent on each other, in that A(m_idx[k][j], k) = m_x[k][j].
	
*/
template <class el_type> 
class skew_symmetry_matrix : public square_matrix<el_type>
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
	skew_symmetry_matrix(int n_rows = 0, int n_cols = 0) : square_matrix<el_type> (n_rows, n_cols)
	{
		sign = -1;
	}
	
	//----Matrix referencing/filling----//
	
	//----Factorizations----//
	
	void ildlpp(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol);
	void ildlmpp(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol);
	void ildlrp(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol);
	
	//------Helpers------//
	
	/*! \brief Updates A.list for iteration k.
		\param k current iteration index.	
	*/
	inline void advance_list(const int& k)
	{
		for (auto it = m_idx[k].begin(); it != m_idx[k].end(); it++)
		{
			ensure_invariant(*it, k, list[*it]); //make sure next element is good.
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
		return list[index].size() + m_idx[index].size();
	}
	
	//----IO helper Functions----//
	// overwrite the helper functions in the class square_matrix

	void put_header(std::string& header)
	{
		header = "%%MatrixMarket matrix coordinate real skew-symmetric";
	}

	inline bool readline(char*& line, int& n_rows, int& n_cols, int& i, int& j, el_type& value)
	{
		sscanf(line, "%d %d %lf", &i, &j, &value);
		i--;
		j--;
		return (i>=0 && j>=0 && i<n_rows&& j<n_cols && i!=j);
	}
};

//------------------ include files for class functions -------------------//

#include "skew_symmetry_matrix_update.h"
#include "skew_symmetry_matrix_ildl.h"
#include "skew_symmetry_matrix_ildlrp.h"
#include "skew_symmetry_matrix_calculate.h"

#endif