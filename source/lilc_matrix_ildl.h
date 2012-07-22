#ifndef _LILC_MATRIX_ILDL_H_
#define _LILC_MATRIX_ILDL_H_


using std::endl;
using std::cout;

template <class el_type>
void lilc_matrix<el_type> :: ildl(lilc_matrix<el_type>& L, block_diag_matrix<el_type>& D, idx_vector_type& perm, double fill_factor, double tol)
{
	int lfil = 2*fill_factor*nnz()/n_cols(); //roughly a factor of 2 since only lower tri. of A is stored
	const double alpha = (1+sqrt(17))/8;  //for use in pivoting.
	el_type w1, wr, d1(0), dr(0), det_D, D_inv11, D_inv22, D_inv12, l_11, l_12;

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
	int i, j, k, r, offset, col_size, col_size2(1);
	bool size_two_piv = false;

	L.resize(ncols, ncols);
	curr_nnzs.reserve(ncols); //makes sure that there is enough space if every element in the column is nonzero
	L.list.resize(ncols ); //allocate a vector of size n for Llist.
	D.resize(ncols );

	for (i = 0; i < ncols; i++) {
		L.m_idx[i].resize(lfil+1);
		L.m_x[i].resize(lfil+1);
	}

	int debug = -1;//16097;
	for (k = 0; k < ncols; k++) {
		size_two_piv = false;

		//zero out work vector
		std::fill (work.begin() + k, work.end(), 0); //just fill the nonzeros. same with temp. otherwise this is O(n^2). totally bad.
		curr_nnzs.clear();

		if (k == debug) { 
			cout << "m_idx[k]: " << m_idx[k] << endl;
			cout << "m_x[k]: " << m_x[k] << endl;
		}
		
		if (m_idx[k].size() > 0) {
			curr_nnzs.assign (m_idx[k].begin(), m_idx[k].end());

			//assigns the non zeros in A(k,:) to the work vector. since only the lower diagonal of A is stored, this is essentially A(k,k+1:n).
			for (j = 0; j < (int) curr_nnzs.size(); j++) {
				work[curr_nnzs[j]] = m_x[k][j];
			}
		}

		//--------------begin pivoting--------------//

		update(k, work, curr_nnzs, L, D, in_set, true);
		
		d1 = work[k];
		work[k] = 0;

		w1 = max(work, curr_nnzs, r);

		if (k == debug) { 
			cout << "curr_nnzs: " << curr_nnzs << endl;
			for (i = 0; i < (int) curr_nnzs.size(); i++) {
				cout << work[curr_nnzs[i]] << " ";
			}
			cout << endl;
			cout << "w1: " << w1 << endl;
		}
		
		if (w1 == 0) {
			//case 0: do nothing. pivot is k.
			if (k == debug) {
				cout << "case 0" << endl;
			}
		} else if (std::abs(d1) >= alpha * w1 ) {
			//case 1: do nothing. pivot is k.
			if (k == debug) {
				cout << "case 1" << endl;
			}
		} else {
			std::fill (temp.begin() + k, temp.end(), 0);
			temp_nnzs.clear();

			ensure_invariant(r, k, list[r], first[r], true);
			offset = first[r];

			for (j = offset; j < (int) list[r].size(); j++) {
				temp_nnzs.push_back(list[r][j]);
				temp[list[r][j]] = coeff(r, list[r][j]);
			}

			unordered_inplace_union(temp_nnzs, m_idx[r].begin(), m_idx[r].end(), in_set);

			for (j = 0; j < (int) m_idx[r].size(); j++) {
				temp[m_idx[r][j]] = m_x[r][j];
			}

			update(r, temp, temp_nnzs, L, D, in_set, true);

			
			dr = temp[r];
			temp[r] = 0;

			wr = max(temp, temp_nnzs, j);

			if (std::abs(d1 * wr)>= alpha*w1*w1) {
				//case 2: do nothing. pivot is k.
				if (k == debug) {
					cout << "case 2" << endl;
					cout << "d1: " << d1 << " wr: " << wr << " w1: " << w1 << endl;
				}
				
				
			} else if (std::abs(dr) >= std::abs(alpha * wr)) {
				//case 3: pivot is k with r: 1x1 pivot case.
				temp[r] = dr;
				work[k] = d1;

				if (k > r) {
					cout << "case 3 " << k << " " << r << endl;
					cout << "fault! " << k+1 << " " << r << endl;
					return;
				}

				//--------pivot A and L ---------//
				pivot(swapk, swapr, swapk_, swapr_, all_swaps, in_set, col_k, col_k_nnzs, col_r, col_r_nnzs, L, k, r);

				//----------pivot rest ----------//
				std::swap(d1, dr);

				//permute perm
				std::swap(perm[k], perm[r]);

				work.swap(temp);	//swap work with temp.
				std::swap(work[k], work[r]);

				curr_nnzs.swap(temp_nnzs);	//swap curr_nnzs with temp_nnzs


				//check if this is working.
				safe_swap(curr_nnzs, k, r);

				//--------end pivot rest---------//

			} else {
				//case 4: pivot is k+1 with r: 2x2 pivot case.
				if (k == debug) {
					cout << "case 4" << " " << k+1 << " " << r << endl;
					cout << "temp_nnzs: " << temp_nnzs << " " << r << endl;
					for (i = 0; i < (int) temp_nnzs.size(); i++) {
						cout << temp[temp_nnzs[i]] << " ";
					}
					cout << endl;
				}
				if (k >= r) {
					cout << "fault! " << k+1 << " " << r << endl;
					return;
				}

				advance_list(k);
				L.advance_first(k);

				temp[r] = dr;
				work[k] = d1;

				size_two_piv = true;

				if (k+1 != r) {
					safe_swap(m_idx[k], k+1, r);
					pivot(swapk, swapr, swapk_, swapr_, all_swaps, in_set, col_k, col_k_nnzs, col_r, col_r_nnzs, L, k+1, r);


					//----------pivot rest ----------//

					//permute perm
					std::swap(perm[k+1], perm[r]);

					//swap two cols of L
					std::swap(work[k+1], work[r]);
					std::swap(temp[k+1], temp[r]);

					safe_swap(curr_nnzs, k+1, r);
					safe_swap(temp_nnzs, k+1, r);
				}
			}
		}
		//--------------end pivoting--------------//

		curr_nnzs.erase(std::remove(curr_nnzs.begin(), curr_nnzs.end(), k), curr_nnzs.end());

		//performs the dual dropping procedure.
		if (!size_two_piv) {
			drop_tol(work, curr_nnzs, lfil, tol);

		} else {
			temp_nnzs.erase(std::remove(temp_nnzs.begin(), temp_nnzs.end(), k), temp_nnzs.end());
			curr_nnzs.erase(std::remove(curr_nnzs.begin(), curr_nnzs.end(), k+1), curr_nnzs.end());
			temp_nnzs.erase(std::remove(temp_nnzs.begin(), temp_nnzs.end(), k+1), temp_nnzs.end());
			
			det_D = d1*dr - work[k+1]*work[k+1];
			//if (det_D != 0) { //assuming matrix is non-singular for now. replace != 0 with EPS later
			D_inv11 = dr/det_D;
			D_inv22 = d1/det_D;
			D_inv12 = -work[k+1]/det_D;
			
			D.off_diagonal(k) = work[k+1];
			D[k+1] = dr;
			
			unordered_inplace_union(curr_nnzs, temp_nnzs.begin(), temp_nnzs.end(), in_set);

			
			for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) { //need -1 on col_size to remove offset from initializing col_size to 1
				l_11 = work[*it]*D_inv11 + temp[*it]*D_inv12;
				l_12 = work[*it]*D_inv12 + temp[*it]*D_inv22;
				
				work[*it] = l_11;
				temp[*it] = l_12;
			}
			
			temp_nnzs.assign(curr_nnzs.begin(), curr_nnzs.end());
			drop_tol(temp, temp_nnzs, lfil, tol);
			drop_tol(work, curr_nnzs, lfil, tol);
			

		}

		D[k] = d1;

		L.m_x[k][0] = 1;
		L.m_idx[k][0] = k;
		count++;

		if (!size_two_piv) {
			i = 0;
			if (k < ncols - 1) //check if D[k] == 0? assuming matrix is non-singular for now.
			for (i = 0; i < (int) curr_nnzs.size(); i++) { //need -1 on col_size to remove offset from initializing col_size to 1
				L.m_idx[k][i+1] = curr_nnzs[i]; //row_idx of L is updated
				L.m_x[k][i+1] = work[curr_nnzs[i]]/D[k]; //work vector is scaled by D[k]

				L.list[curr_nnzs[i]].push_back(k); //update Llist
				count++;
			}
			
			col_size = 1 + i;
			//update lfirst
			L.advance_first(k);
			advance_list(k);
		} else {
			
			D.off_diagonal(k) = work[k+1];
			D[k+1] = dr;

			L.m_x[k+1][0] = 1;
			L.m_idx[k+1][0] = k+1;
			count++;

			i = 0;
			for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
				if (work[*it] != 0) {
					L.m_x[k][i+1] = work[*it]; //row_idx of L is updated
					L.m_idx[k][i+1] = *it; //work vector is scaled by D[k]
					L.list[*it].push_back(k); //update Llist
					count++;
					i++;
				}
			}
			
			j = 0;
			for (auto it = temp_nnzs.begin(); it != temp_nnzs.end(); it++) {
				if (temp[*it] != 0) {
					L.m_x[k+1][j+1] = temp[*it]; //row_idx of L is updated
					L.m_idx[k+1][j+1] = *it; //work vector is scaled by D[k]
					L.list[*it].push_back(k+1); //update Llist
					count++;
					j++;
				}
				
			}

			col_size = 1 + i;
			col_size2 = 1 + j;

			//update lfirst
			L.advance_first(k+1);
			advance_list(k+1);
			
			temp[k] = 0;
			temp[k+1] = 0;
			work[k+1] = 0;
			for (auto it = temp_nnzs.begin(); it != temp_nnzs.end(); it++) {
				temp[*it] = 0;
			}
			//}
		}
		
		//cleanup
		work[k] = 0;
		for (auto it = curr_nnzs.begin(); it != curr_nnzs.end(); it++) {
			work[*it] = 0;
		}
		
		L.m_x[k].resize(col_size);
		L.m_idx[k].resize(col_size);

		if (k == debug) {
			cout << "part of two piv? " << size_two_piv << endl;
			cout << "lfil: " << lfil << endl;
			cout << "D[k]: " << D[k] << endl;
			cout << "D(k,k+1): " << D.off_diagonal(k) << endl;
			cout << "D[k+1]: " << D[k+1] << endl;
			
			cout << "curr_nnzs: " << curr_nnzs << " " << r << endl;
			for (i = 0; i < col_size; i++) {
				cout << work[curr_nnzs[i]]/D[k] << " ";
			}
			cout << endl;
		}
		
		if (size_two_piv) {
			L.m_x[k+1].resize(col_size2);
			L.m_idx[k+1].resize(col_size2);
			k++;
		}
		
		// std::cout << k << std::endl;
		// std::string s;
		// std::cin >> s;
		// if (s == "l") {
		// std::cout << L << std::endl;
		// } else if (s == "a") {
		// std::cout << to_string() << std::endl;
		// }
		
	}

	L.nnz_count = count;

}

#endif
