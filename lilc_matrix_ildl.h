#ifndef _LILC_MATRIX_ILDL_H_
#define _LILC_MATRIX_ILDL_H_

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
	elt_vector_type work(ncols, 0), temp(ncols, 0), col_k, col_r;
	idx_vector_type curr_nnzs, temp_nnzs, col_k_nnzs, col_r_nnzs, all_swaps;  //non-zeros on current col.
	vector<bool> in_set(ncols, 0);
	std::pair<idx_it, elt_it> its_k, its_r;
	vector<idx_it> swapk, swapr;
	vector<list_it> swapk_, swapr_;
	
	int count = 0; //the current non-zero in L.
	int i, j, k, r, offset, col_size, col_size2;
	bool size_two_piv = false;

	L.resize(ncols, ncols);
	curr_nnzs.reserve(ncols); //makes sure that there is enough space if every element in the column is nonzero
	L.list.resize(ncols ); //allocate a vector of size n for Llist.
	D.resize(ncols ); 
	
	for (i = 0; i < ncols; i++) {
		L.m_idx[i].resize(lfil+1);
		L.m_x[i].resize(lfil+1);
	}
	
	for (k = 0; k < ncols; k++) {
		size_two_piv = false;
		
		//zero out work vector
		std::fill (work.begin() + k, work.end(), 0);
		col_size = 1;
		
		//future self: remember you need m_idx[k]+1 only if there is a diag elem in A.col(k)
		if (m_idx[k].size() > 0) {
			offset = (m_idx[k][0] == k ? 1 : 0);
			curr_nnzs.assign (m_idx[k].begin()+offset, m_idx[k].end());
			
			//assigns the non zeros in A(k,:) to the work vector. since only the lower diagonal of A is stored, this is essentially A(k,k+1:n).
			for (j = 0; j < (int) curr_nnzs.size(); j++) {
				work[curr_nnzs[j]] = m_x[k][j+offset];
			}
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
					
				} else if (std::abs(dr) >= std::abs(alpha * wr)) {
					//case 3: pivot is k with r: 1x1 pivot case.
					//--------pivot A and L ---------//
					pivot(swapk, swapr, swapk_, swapr_, all_swaps, in_set, col_k, col_k_nnzs, col_r, col_r_nnzs, L, k, r);

					//----------pivot rest ----------//
					temp[r] = dr;
					work[k] = d1;
					std::swap(d1, dr);
					
					//permute perm
					std::swap(perm[k], perm[r]);
					
					work.swap(temp);	//swap work with temp.
					std::swap(work[k], work[r]);
					
					curr_nnzs.swap(temp_nnzs);	//swap curr_nnzs with temp_nnzs
					
					for (auto it = curr_nnzs.begin(); it!= curr_nnzs.end(); it++) {
						if (*it == r && work[r] == 0) 
						curr_nnzs.erase(it);
					}
					//--------end pivot rest---------//
					
				} else {
					continue;
					//case 4: pivot is k+1 with r: 2x2 pivot case.
					size_two_piv = true;
					
					pivot(swapk, swapr, swapk_, swapr_, all_swaps, in_set, col_k, col_k_nnzs, col_r, col_r_nnzs, L, k+1, r);

					//----------pivot rest ----------//
					temp[r] = dr;
					work[k] = d1;
					std::swap(d1, dr);
					
					//permute perm
					std::swap(perm[k+1], perm[r]);
					
					//swap two cols of L
					std::swap(work[k+1], work[r]);
					std::swap(temp[k+1], temp[r]);
					
					safe_swap(curr_nnzs, k+1, r);
					safe_swap(temp_nnzs, k+1, r);
				}
			}
			//--------------end pivoting--------------//
			//update lfirst
			for (auto it = L.list[k].begin(); it != L.list[k].end(); it++) {
				L.first[*it]++;
			}
			
			if (m_idx[k].size() > 0) {
				offset = (m_idx[k][0] == k ? 1 : 0);
				for (j = offset; j < (int) m_idx[k].size(); j++) {
					if (!list[m_idx[k][j]].empty())
					list[m_idx[k][j]].pop_front();
				}
			}
			
			if (size_two_piv) {
				if (m_idx[k+1].size() > 0) {
					offset = (m_idx[k+1][0] == k+1 ? 1 : 0);
					for (j = offset; j < (int) m_idx[k+1].size(); j++) {
						if (!list[m_idx[k+1][j]].empty())
						list[m_idx[k+1][j]].pop_front();
					}
				}
			}
			//performs the dual dropping procedure.
			drop_tol(work, curr_nnzs, lfil, tol);
			col_size += std::min(lfil, (int) curr_nnzs.size());
			
			if (size_two_piv) {
				drop_tol(temp, temp_nnzs, lfil, tol);
				col_size2 = 1+std::min(lfil, (int) temp_nnzs.size());
			}
		}
		
		if (!size_two_piv) {
			D[k] = d1;
			
			L.m_x[k][0] = 1;
			L.m_idx[k][0] = k;
			count++;
			
			if (D[k] == 0) col_size = 1;
			if (k < ncols - 1)
				for (i = 0; i < col_size-1; i++) //need -1 on col_size to remove offset from initializing col_size to 1
				{
					L.m_idx[k][i+1] = curr_nnzs[i]; //row_idx of L is updated
					L.m_x[k][i+1] = work[curr_nnzs[i]]/D[k]; //work vector is scaled by D[k]
					
					L.list[curr_nnzs[i]].push_back(k); //update Llist
					count++;
				}
		}
		
		
		L.m_x[k].resize(col_size);
		L.m_idx[k].resize(col_size);
	}
	
	L.nnz_count = count;
	
}

#endif 
