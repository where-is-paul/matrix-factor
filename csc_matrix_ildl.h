#ifndef _CSC_MATRIX_ILDL_H_
#define _CSC_MATRIX_ILDL_H_

#include <algorithm>
#include <cmath>
#include <list>
#include <deque>
#include "csc_matrix_ildl_helpers.h"

/*! possible optimizations/to-do's: 
	- some easy parallelization of for loops. 
	
	- using a vector to accumulate all nonzeros in previous cols
	  during iteration.
	  
	- small changes to loops (e.g. assign .begin(), and .end() so that they will not be evaluated each iter.
	---> for (auto itr = new_values.begin(), end_itr = new_values.end(); itr != end_itr; ++itr ) )
	
	- make a special diagonal matrix class for storing 1x1 and 2x2 pivots later.
	
	- typedef some of the std::vector<...> into more readable names.
	
	- use epsilon tolerances in comparing doubles later.
	
*/

template <class idx_type, class el_type>
void csc_matrix<idx_type, el_type> :: ildl(csc_matrix<idx_type, el_type>& L, elt_vector_type& D, int lfil, double tol)
{	
	//Llist is a deque of linked lists that gives the non-zero elements in each row of L. since at any time we may swap between two rows, we require a linked lists for each row of L. a deque is used as it might be desireable to deallocate all linked lists for rows i < k on step k (this is currently not done, as the memory used in maintaining linked lists for all rows is not much).
	std::vector< std::deque< idx_type > > Llist; 
	
	//work is a work vector for the current column. Lfirst is a linked list that gives the first nonzero element in column k with row index i > k. (i.e. the first nonzero in L(k+1:n, k).
	std::vector<el_type> work(n_cols(), 0), temp(n_cols(), 0);
	
	std::vector<idx_type> curr_nnzs, temp_nnzs, Lfirst(n_cols(), 0); //non-zeros on current col.
	
	int count = 0; //the current non-zero in L.
	const double alpha = (1+sqrt(17))/8;  //for use in pivoting.
	idx_type i, j, k, offset, r;
	el_type l_ik, w1, wr;
	
	curr_nnzs.reserve(n_cols()); //makes sure that there is enough space if every element in the column is nonzero
	Llist.resize(n_cols()); //allocate a vector of size n for Llist.
	
	L.resize(n_rows(), n_cols(), (lfil+1)*n_cols()); //(+1 because there are 1s on the diagonal. they wont need to be stored if we want to optimize)
	D.resize(n_cols()); 
	
	for (i = 0; i < n_cols(); i++) {
		D[i] = coeff(i,i);
	}
	
	for (k = 0; k < n_cols(); k++) {
	    //zero out work vector
	    std::fill (work.begin() + k, work.end(), 0);
	    
	    //the +1 avoids assigning diagonal element as nonzero since its stored in D. the min is in case m_col_idx[k] == m_col_idx[k+1] (an empty column).
	    curr_nnzs.assign (m_row_idx.begin() + min(m_col_idx[k] + 1, m_col_idx[k+1]), m_row_idx.begin() + m_col_idx[k+1]);
		
		//assigns the non zeros in A(k,:) to the work vector. since only the lower diagonal of A is stored, this is essentially A(k,k+1:n).
	    for (j = m_col_idx[k] + 1; j < m_col_idx[k+1]; j++) {
			work[m_row_idx[j]] = m_x[j];
		}
		
		if (k < n_cols() - 1) {
			//--------------begin pivoting--------------//
			//perform delayed updates on the current col (k) of A
			update(k, work, curr_nnzs, L, D, Lfirst, Llist);
			
			w1 = max(work, curr_nnzs, r);
			if (w1 == 0) {
				//case 0: do nothing. pivot is k.
			} else if (std::abs(D[k]) >= alpha * w1 ) {
				//case 1: do nothing. pivot is k.
			} else {
				continue; //have not finished implementing pivoting just yet.
				
				//assign A(k+1:n, r) to temp. also update the diagonal elements. we'll need a_rr
				/* 
				
					fill in code here 
				
				*/
				
				//perform delated updates on col (r) of A.
				update(r, temp, temp_nnzs, L, D, Lfirst, Llist);
				wr = max(temp, temp_nnzs, r);
				if (std::abs(D[k] * wr)>= alpha*w1*w1) {
					//case 2: do nothing. pivot is k.
				} else if (std::abs(D[r]) >= alpha * wr) {
					//case 3: pivot is k with r. 1x1 pivot.
					/*
					pivot(k,r);
					*/
					
					work.swap(temp);	//swap work with temp.
					curr_nnzs.swap(temp_nnzs);	//swap curr_nnzs with temp_nnzs
				} else {
					//case 4: pivot is k+1 with r: 2x2 pivot case.
					/*
					pivot(k+1,r);
					*/
				}
			}
			
			//--------------end pivoting--------------//
			
			//performs the dual dropping procedure.
			drop_tol(work, curr_nnzs, lfil, tol);
			
		}
		//get 1s on the diagonal
		L.m_row_idx[count] = k;
		L.m_x[count] = 1;
		count++;

		if (k < n_cols() - 1)
		for (i = 0; i < (idx_type) std::min(lfil, (int) curr_nnzs.size()); i++) {
		    L.m_row_idx[count] = curr_nnzs[i]; //row_idx of L is updated
		    L.m_x[count] = work[curr_nnzs[i]]/D[k]; //work vector is scaled by D[k]
			Llist[curr_nnzs[i]].push_back(k); //update Llist
			count++;
		}
		
		L.m_col_idx[k+1] = count; //the end of the current column is assigned to col_idx
		
		if (k < n_cols() - 1) {
			//finds out where L(k+1:n, k) starts
			offset = L.m_col_idx[k] + Lfirst[k] + 1;
			for (i = offset; i < L.m_col_idx[k+1]; i++) {
				l_ik = L.m_x[i];
				D[L.m_row_idx[i]] -= l_ik * D[k] * l_ik;	//update diagonal
			}
		}
	}
	
	//resize vectors of L down to the right size.
	L.m_row_idx.resize(count); 
	L.m_x.resize(count);
	
}

#endif 
