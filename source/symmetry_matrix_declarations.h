// -*- mode: c++ -*-
#ifndef SYMMETRY_MATRIX_DECLARATIONS_H
#define SYMMETRY_MATRIX_DECLARATIONS_H

#include "symmetry_swap_struct.h"
#include "ultriangular.h"

/*! \brief A list-of-lists (LIL) matrix in column oriented format.

	For convience, the matrix this class represents will be refered to as matrix A.
	In LIL-C format, each column of A (an n*n matrix) is stored as a separate vector. The nonzeros are stored in m_idx while the non-zeros are stored in m_x. Both m_x and m_idx are initialized to a list of n lists. m_idx and m_x are ordered dependent on each other, in that A(m_idx[k][j], k) = m_x[k][j].
	
*/
template <class el_type> 
class symmetry_matrix : public half_matrix<el_type>
{
public:
	
	//-------------- typedefs and inherited variables --------------//
	using half_matrix<el_type>::m_idx; // linux compiler gcc need these
	using half_matrix<el_type>::m_x;
	using half_matrix<el_type>::m_n_rows;
	using half_matrix<el_type>::m_n_cols;
	using half_matrix<el_type>::n_rows;
	using half_matrix<el_type>::n_cols;
	using half_matrix<el_type>::nnz;
	using half_matrix<el_type>::nnz_count;

	// using half_matrix<el_type>::idx_vector_type;
	// using half_matrix<el_type>::elt_vector_type;
	// using half_matrix<el_type>::idx_it;
	// using half_matrix<el_type>::elt_it;

	using half_matrix<el_type>::list;
	using half_matrix<el_type>::list_first;

	using half_matrix<el_type>::S;

	typedef typename half_matrix<el_type>::idx_vector_type idx_vector_type; // linux compiler gcc need these
	typedef typename half_matrix<el_type>::elt_vector_type elt_vector_type;
	typedef typename half_matrix<el_type>::idx_it idx_it;
	typedef typename half_matrix<el_type>::elt_it elt_it;

	// typedef typename idx_vector_type::iterator idx_it;
	// typedef typename elt_vector_type::iterator elt_it;
	
public:
	
	/*! \brief Constructor for a column oriented list-of-lists (LIL) matrix. Space for both the values list and the indices list of the matrix is allocated here.
	*/
	symmetry_matrix(int n_rows = 0, int n_cols = 0) : half_matrix<el_type> (n_rows, n_cols)
	{
	}
	
	//----Matrix referencing/filling----//
	
	/*! \brief Finds the (i,j)th coefficient of the matrix.
		\param i the row of the (i,j)th element (zero-indexed).
		\param j the col of the (i,j)th element (zero-indexed).
		\param offset an optional search offset for use in linear search (start at offset instead of 0).
		\return The (i,j)th element of the matrix.
	*/
	inline el_type coeff(const int& i, const int& j, int offset = 0) const 
	{	
		for (unsigned int k = offset, end = m_idx[j].size(); k < end; k++)
			if (m_idx[j][k] == i)
				return m_x[j][k];
		
		return 0;
	}
	
	//-----Reorderings/Rescalings------//


	/*! \brief Given a permutation vector perm, A is permuted to P'AP, where P is the permutation matrix associated with perm. 
		\param perm the permutation vector.
	*/
	void sym_perm(vector<int>& perm);
	
	/*!	\brief The symmetric matrix A is equilibrated and the symmetric equilibrated matrix SAS is stored in A, where S is a diagonal scaling matrix. 
		
		This algorithm is based on the one outlined in "Equilibration of Symmetric Matrices in the Max-Norm" by Bunch (1971).
	*/
	void sym_equil();
	
	//----Factorizations----//
	/*! \brief Performs an LDL' factorization of this matrix. 
		
		The pivoted matrix P'AP will be stored in place of A. In addition, the L and D factors of P'AP will be stored in L and D (so that P'AP = LDL'). The factorization is performed in crout order and follows the algorithm outlined in "Crout versions of the ILU factorization with pivoting for sparse symmetric matrices" by Li and Saad (2005).
	
		\param L the L factor of this matrix.
		\param D the D factor of this matrix.
		\param perm the current permutation of A.
		\param fill_factor a parameter to control memory usage. Each column is guaranteed to have fewer than fill_factor*(nnz(A)/n_col(A)) elements.
		\param tol a parameter to control agressiveness of dropping. In each column, elements less than tol*norm(column) are dropped.
	*/
	void ildl(ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, const double& fill_factor, const double& tol);
	
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
	/*! \brief Performs a symmetric permutation between row/col k & r of A.
	
		\param s a struct containing temporary variables needed during pivoting.
		\param in_set a bitset needed for unordered unions during pivoting.
		\param L the lower triangular factor of A.
		\param k index of row/col k.
		\param r index of row/col r.
	*/
	inline void pivot(symmetry_swap_struct<el_type> s, ultriangular_matrix<el_type>& L, const int& k, const int& r);
	
	/*! \brief Ensures two the invariants observed by A.first and A.list are held.
		
		\invariant
		If this matrix is a lower triangular factor of another matrix:
			-# On iteration k, first[i] will give the number of non-zero elements on col i of A before A(k, i).
			-# On iteration k, list[i][ first[i] ] will contain the first element below or on index k of column i of A.
		
		\invariant
		If this matrix is the matrix to be factored:
			-# On iteration k, first[i] will give the number of non-zero elements on row i of A before A(i, k).
			-# On iteration k, list[i][ first[i] ] will contain the first element right of or on index k of row i of A.
			
		\param j the column of con.
		\param k the iteration number.
		\param con the container to be swapped.
		\param update_list boolean indicating whether list or m_x/m_idx should be updated.
	*/
	inline void ensure_invariant(const int& j, const int& k, idx_vector_type& con)
	{
		int offset = list_first[j];
		if (con[offset] == k)
			return;
		
		int i, min = offset;
		for (i = offset; i < (int) con.size(); i++)
		{
			if (con[i] == k) {
				min = i; 
				break;
			}
		}
		
		std::swap(con[offset], con[min]);
	}
	
	/*! \brief Updates A.list for iteration k.
		\param k current iteration index.	
	*/
	inline void advance_list(const int& k)
	{
		for (auto it = m_idx[k].begin(); it != m_idx[k].end(); it++)
		{
			if (*it == k) continue;
			ensure_invariant(*it, k, list[*it]); //make sure next element is good.
			list_first[*it]++; //invariant ensured.
		}
	}
	
	inline void update_single(const int& j, const el_type& l_ki, const el_type& d, std::vector<el_type>& work, std::vector<int>& curr_nnzs, ultriangular_matrix<el_type>& L, vector<bool>& in_set);
	inline void update(const int& r, std::vector<el_type>& work, std::vector<int>& curr_nnzs, ultriangular_matrix<el_type>& L, block_diag_matrix<el_type>& D, vector<bool>& in_set);

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
	// overwrite the helper functions in the class half_matrix

	void put_header(std::string& header)
	{
		header = "%%MatrixMarket matrix coordinate real symmetric";
	}

	inline bool readline (char*& line, int& n_rows, int& n_cols, int& i, int& j, el_type& value)
	{
		sscanf(line, "%d %d %lf", &i, &j, &value);
		i--;
		j--;
		return (i>=0 && j>=0 && i<n_rows&& j<n_cols);
	}
};

//------------------ include files for class functions -------------------//

#include "symmetry_matrix_perm.h"
#include "symmetry_matrix_equil.h"
#include "symmetry_matrix_ildl.h"
#include "symmetry_matrix_ildlrp.h"
#include "symmetry_matrix_update.h"
#include "symmetry_matrix_pivot.h"

#endif