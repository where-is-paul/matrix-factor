#ifndef _LILC_MATRIX_ILDL_H_
#define _LILC_MATRIX_ILDL_H_

#include <algorithm>
#include <cmath>
#include <list>
#include <deque>
#include "lilc_matrix_ildl_helpers.h"

template <class el_type>
void lilc_matrix<el_type> :: ildl(lilc_matrix<el_type>& L, elt_vector_type& D, idx_vector_type& perm, int lfil, double tol)
{	
	const double alpha = (1+sqrt(17))/8;  //for use in pivoting.
	el_type w1, wr, d1, dr;
	
	//L.list is a deque of linked lists that gives the non-zero elements in each row of L. 
	//since at any time we may swap between two rows, we require a linked lists for each row of L. 
	//A deque is used as it might be desirable to deallocate all linked lists for rows i < k on step k 
	//(this is currently not done, as the memory used in maintaining linked lists for all rows is not much).
	
	const int ncols = n_cols();
	//work is a work vector for the current column. L.first is a linked list that gives the first nonzero element in column k with row index i > k. (i.e. the first nonzero in L(k+1:n, k).
	elt_vector_type work(ncols , 0), temp(ncols , 0), col_k, col_r;
	idx_vector_type curr_nnzs, temp_nnzs, col_k_nnzs, col_r_nnzs, all_swaps;  //non-zeros on current col.
	vector<bool> in_set(ncols , 0);
	std::pair<idx_it, elt_it> its_k, its_r;
	vector<idx_it> swapk, swapr;
	vector<list_it> swapk_, swapr_;
	
	int count = 0; //the current non-zero in L.
	int i, j, k, r, offset, col_size;

	L.resize(ncols ,ncols);
	curr_nnzs.reserve(ncols ); //makes sure that there is enough space if every element in the column is nonzero
	L.list.resize(ncols ); //allocate a vector of size n for Llist.
	D.resize(ncols ); 
	
	for (i = 0; i < ncols; i++) {
		L.m_idx[i].resize(lfil+1);
		L.m_x[i].resize(lfil+1);
	}
	
	for (k = 0; k < ncols; k++) {
				
	    //zero out work vector
	    std::fill (work.begin() + k, work.end(), 0);
		col_size = 1;
	    
	    //future self: remember you need m_idx[k]+1 only if there is a diag elem in A.col(k)
	    offset = (m_idx[k][0] == k ? 1 : 0);
	    curr_nnzs.assign (m_idx[k].begin()+offset, m_idx[k].end());
		
		//assigns the non zeros in A(k,:) to the work vector. since only the lower diagonal of A is stored, this is essentially A(k,k+1:n).
	    for (j = 0; j < (int) curr_nnzs.size(); j++) {
			work[curr_nnzs[j]] = m_x[k][j+offset];
		}
		
		d1 = coeff(k,k);
		for (auto it = L.list[k].begin(); it != L.list[k].end(); it++) { 
			offset = L.first[*it];
			d1 -= L.m_x[*it][offset] * D[*it] * L.m_x[*it][offset];	//update diagonal
		}
		
		if (k < ncols - 1) {
			//--------------begin pivoting--------------//
			update(k, work, curr_nnzs, L, D, in_set);
			
			w1 = max(work, curr_nnzs, r);
			if (w1 == 0) {
				//case 0: do nothing. pivot is k.
			} else if (std::abs(d1) >= alpha * w1 ) {
				//case 1: do nothing. pivot is k.
			} else {
				std::fill (temp.begin() + k, temp.end(), 0);
				temp_nnzs.clear();
				
				temp[k] = coeff(k,r);
				for (auto it = L.list[k].begin(); it != L.list[k].end(); it++) { 
					offset = L.first[*it];
					temp[k] -= L.coeff(r, *it) * D[*it] * L.m_x[*it][offset];
				}
				
				offset = (list[r][0] == k ? 1 : 0);
				for (j = 0; j < (int) list[r].size(); j++) {
					temp_nnzs.push_back(list[r][j]);
					temp[list[r][j]] = coeff(r, list[r][j]);
				}
				
				unordered_inplace_union(temp_nnzs, m_idx[r].begin(), m_idx[r].end(), in_set);
				for (j = 0; j < (int) temp_nnzs.size(); j++) 
					in_set[temp_nnzs[j]] = false;
				
				for (j = 0; j < (int) m_idx[r].size(); j++) {
					temp[m_idx[r][j]] = m_x[r][j];
				}
				
				update(r, temp, temp_nnzs, L, D, in_set, true);
				
				dr = temp[r];
				temp[r] = 0;
				
				wr = max(temp, temp_nnzs, j);
				
				if (std::abs(d1 * wr)>= alpha*w1*w1) {
					//case 2: do nothing. pivot is k.
					
				} else if (true && std::abs(dr) >= std::abs(alpha * wr)) {
					//case 3: pivot is k with r: 1x1 pivot case.
					col_k.clear();
					col_k_nnzs.clear();
					
					col_r.clear();
					col_r_nnzs.clear();
					
					temp[r] = dr;
					work[k] = d1;
					
					//----------pivot A ----------//
					swapr_.clear();
					swapk_.clear();
					all_swaps.clear();
					
					//swap A(:, k) with A(:, r)
					
					//invariant: ensure m_idx[:][0] and m_x[:][0] contains the smallest index element.
					if (coeff(r, r) !=0){
						col_k_nnzs.push_back(k);
						col_k.push_back(coeff(r, r));
					}
					
					if (coeff(k, k) !=0){
						col_r_nnzs.push_back(r);
						col_r.push_back(coeff(k, k));
					}
					
					for (i = 0; i < (int) list[r].size(); i++) {
						coeffRef(r, list[r][i], its_k);
						col_k_nnzs.push_back(*its_k.first);
						col_k.push_back(*its_k.second);
						
						*its_k.first = m_idx[list[r][i]].back();
						*its_k.second = m_x[list[r][i]].back();
						
						m_idx[list[r][i]].pop_back();
						m_x[list[r][i]].pop_back();
					}
					
					offset = (m_idx[r][0] == r ? 1 : 0);
					std::copy(m_x[r].begin()+offset, m_x[r].end(), std::back_inserter(col_k));
					std::copy(m_idx[r].begin()+offset, m_idx[r].end(), std::back_inserter(col_k_nnzs));
					
					//invariant: ensure list[:][0] contain the index nearest in value to k.
					for (auto it = m_idx[r].begin() + offset; it != m_idx[r].end(); it++) {
						for (i = 0; i < (int) list[*it].size(); i++) {
							if (list[*it][i] == r) {
								swapk_.push_back(list[*it].begin() + i);
								all_swaps.push_back(*it);
								break;
							}
						}
					}	
					
					//swap A(k:r, k) with A(r, k:r);
					offset = (m_idx[k][0] == k ? 1 : 0);
					for (i = 0; i < (int) m_idx[k].size(); i++) {
						if (m_idx[k][i] < r) {
							m_idx[m_idx[k][i]].push_back(r);
							m_x[m_idx[k][i]].push_back(m_x[k][i]);
						} else if (m_idx[k][i] != r) {
							col_r.push_back(m_x[k][i]);
							col_r_nnzs.push_back(m_idx[k][i]);
							
							//this part can be simplified since invariant ensures list[:][0] is always the idx closest to k. wont make a diff on running time though
							for (j = 0; j < (int) list[m_idx[k][i]].size(); j++) {
								if (list[m_idx[k][i]][j] == k) {
									swapr_.push_back(list[m_idx[k][i]].begin() + j);
									all_swaps.push_back(m_idx[k][i]);
									break;
								}
							}
						}
					}
					
					for (auto it = swapk_.begin(); it != swapk_.end(); it++) {
						**it = k;
					}
					
					for (auto it = swapr_.begin(); it != swapr_.end(); it++) {
						**it = r;
					}
					
					int min = 0;
					for (auto it = all_swaps.begin(); it != all_swaps.end(); it++) {
						for (i = 0; i < (int) list[*it].size(); i++) {
							min = 0;
							if (list[*it][i] == k) {
								min = i; break;
							} else if ( list[*it][i] < list[*it][min] ) {
								min = i;
							}
						}
						
						std::swap(list[*it][0], list[*it][min]);
						std::swap(list[*it][0], list[*it][min]);
					}
					
					//set the kth col
					m_idx[k].swap(col_k_nnzs);
					m_x[k].swap(col_k);
					
					//set the rth col
					m_idx[r].swap(col_r_nnzs);
					m_x[r].swap(col_r);
					
					list[k].swap(list[r]);
					//--------end pivot A---------//
					
					//----------pivot L ----------//
					swapr.clear();
					swapk.clear();
					all_swaps.clear();
					
					for (auto it = L.list[k].begin(); it != L.list[k].end(); it++)
					{
						for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
							if (L.m_idx[*it][i] == k) {
								swapr.push_back(L.m_idx[*it].begin() + i);
								break;
							}
						}
					}
					
					for (auto it = L.list[r].begin(); it != L.list[r].end(); it++) {
						for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
							if (L.m_idx[*it][i] == r) {				
								swapk.push_back(L.m_idx[*it].begin() + i);
								break;
							}
						}
					}
					
					for (auto it = swapk.begin(); it != swapk.end(); it++) {
						**it = k;
					}
					
					for (auto it = swapr.begin(); it != swapr.end(); it++) {
						**it = r;
					}
					
					all_swaps.assign(L.list[r].begin(), L.list[r].end());
					unordered_inplace_union(all_swaps, L.list[k].begin(), L.list[k].end(), in_set);
					
					//invariant: ensure m_idx[:][L.first[i]+1] all contain the index nearest in value to k. 
					min = 0;
					for (auto it = all_swaps.begin(); it != all_swaps.end(); it++) {
						in_set[*it] = 0; //reset in_set for future use;
						for (i = L.first[*it]; i < (int) L.m_idx[*it].size(); i++) {
							min = L.first[*it];
							if (L.m_idx[*it][i] == k) {
								min = i; break;
							} else if ( L.m_idx[*it][i] < L.m_idx[*it][min] ) {
								min = i;
							}
						}
						
						std::swap(L.m_idx[*it][L.first[*it]], L.m_idx[*it][min]);
						std::swap(L.m_x[*it][L.first[*it]], L.m_x[*it][min]);
					}
					
					L.list[k].swap(L.list[r]);
					//--------end pivot L---------//
					
					//----------pivot D ----------//
					std::swap(D[k], D[r]);
					std::swap(d1, dr);
					//--------end pivot D---------//
					
					//----------pivot rest ----------//
					//permute perm
					std::swap(perm[k], perm[r]);
					
					//swap k and r in work later
					work.swap(temp);	//swap work with temp.
					std::swap(work[k], work[r]);
					
					curr_nnzs.swap(temp_nnzs);	//swap curr_nnzs with temp_nnzs
					
					for (auto it = curr_nnzs.begin(); it!= curr_nnzs.end(); it++) {
						if (*it == r && work[r] == 0) 
							curr_nnzs.erase(it);
					}
					//--------end pivot rest---------//
					
				} else {
					//case 4: pivot is k+1 with r: 2x2 pivot case.
					
					
				}
			}
			
			//--------------end pivoting--------------//
			
			//update lfirst
			for (auto it = L.list[k].begin(); it != L.list[k].end(); it++) {
				L.first[*it]++;
			}
		
			offset = (m_idx[k][0] == k ? 1 : 0);
			for (j = offset; j < (int) m_idx[k].size(); j++) {
				if (!list[m_idx[k][j]].empty())
					list[m_idx[k][j]].pop_front();
				//first[m_idx[k][j]]++;
			}
			
			//performs the dual dropping procedure.
			drop_tol(work, curr_nnzs, lfil, tol);
			col_size += std::min(lfil, (int) curr_nnzs.size());
			
		}
		
		D[k] = d1;
				
		L.m_x[k][0] = 1;
		L.m_idx[k][0] = k;
		count++;
		
		if (k < ncols - 1 && D[k] != 0)
		for (i = 0; i < col_size-1; i++) //need -1 on col_size to remove offset from initializing col_size to 1
		{
		    L.m_idx[k][i+1] = curr_nnzs[i]; //row_idx of L is updated
		    L.m_x[k][i+1] = work[curr_nnzs[i]]/D[k]; //work vector is scaled by D[k]
		    
			L.list[curr_nnzs[i]].push_back(k); //update Llist
			count++;
		}
		
		L.m_x[k].resize(col_size);
		L.m_idx[k].resize(col_size);
		
	}
	
	L.nnz_count = count;
	
}

#endif 
