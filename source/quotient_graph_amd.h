
//overall todo list:
//	change variable names to make more sense.
//  documentation
//  remove explicit using namespace, etc.
#include <queue>
using namespace std;

#ifndef DEBUG
#define DEBUG
template<class Container>
std::ostream& operator<< (std::ostream& os, const Container& vec)
{
	os << "[";
	if (!vec.empty())
	{
		for (auto it = vec.begin(); it+1 != vec.end(); it++)
		{
			os << *it << ", ";
		}
		
		it++;
		os << *it;
	}
	os << "]";
	return os;
}
#endif

class quotient_graph_amd {
	typedef vector<vector<int> > adj_list;
	
	adj_list graph;
	vector<bool> in_set, found, eliminated;
	vector<int> reach, quo_nodes, temp, degree, curr_nbrs, temp_nbrs, L_approx;
	int num_elim;
	
	public:
		quotient_graph_amd(adj_list& g) {
			int n = g.size();
			in_set.resize(n, false);
			found.resize(n, false);
			eliminated.resize(n, false);
			degree.resize(n);
			L_approx.resize(n, -1);
			
			temp.reserve(n); 
			reach.reserve(n);
			quo_nodes.reserve(n);
			graph.resize(n);
			curr_nbrs.reserve(n);
			temp_nbrs.reserve(n);
			
			num_elim = 0;
			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < (int) g[i].size(); j++) {
					int v = g[i][j];
					if (v == i) continue;
					graph[i].push_back(v);
					graph[v].push_back(i);
				}
			}
			
			for (int i = 0; i < n; i++) {
				degree[i] = graph[i].size();
			}
		}
		
		inline int min_deg() {
			int best = 0, deg_best = 10*degree.size();
			for (int i = 0; i < (int) degree.size(); i++) {
				if (eliminated[i]) continue;
				if (degree[i] < deg_best) {
					best = i;
					deg_best = degree[i];
				}
			}
			
			return best;
		}
		
		inline void update_temp(const int& i) {
			temp.clear();
			
			int n = graph[i].size(), m = quo_nodes.size();
			for (int j = 0; j < n; j++) {
				in_set[graph[i][j]] = true;
			}
			
			for (int j = 0; j < m; j++) {
				if (in_set[quo_nodes[j]]) 
					temp.push_back(quo_nodes[j]);
			}
			
			for (int j = 0; j < n; j++) {
				in_set[graph[i][j]] = false;
			}
		}
		
		inline void reachable(const int& i) {
			reach.clear();
			
			int w, k; queue<int> q;
			for (int j = 0; j < (int) graph[i].size(); j++) {
				w = graph[i][j];
				if (eliminated[w]) {
					q.push(w);
				} 
				in_set[w] = true;
			}
		
			in_set[i] = true;
			found[i] = true;
			while (!q.empty()) {
				k = q.front(); q.pop();

				found[k] = true;
				for (int j = 0; j < (int) graph[k].size(); j++) {
					w = graph[k][j];
					if (w == i) {
						swap(graph[k][j], graph[k][graph[k].size()-1]);
						graph[k].pop_back();
						if (j < (int) graph[k].size())
							w = graph[k][j];
					}
					
					//the _found_ bitset serves dual purposes here:
					//it prevents the bfs from traversing repeated nodes
					//while also prevents the loop below from adding
					//repeated notes to reach.
					while (j < (int) graph[k].size() && in_set[w] && !found[w]) {
						found[w] = true;
						reach.push_back(w);
						swap(graph[k][j], graph[k][graph[k].size()-1]);
						graph[k].pop_back();
						if (j < (int) graph[k].size())
							w = graph[k][j];
					}

					if (j >= (int) graph[k].size()) break;
					
					if (eliminated[w]) {
						if (!found[w])
							q.push(w);
					} else if (!in_set[w]) {
						reach.push_back(w);
						found[w] = true;
					}
					
					//the in_set bitset is a marker for removing repeated
					//nodes from the adj. list of the supernodes.
					in_set[w] = true;
				}
			}		
			
			//cleaning up...
			//the two lines below can be optimized by keeping track of exactly what was put in.
			fill(in_set.begin(), in_set.end(), false);
			fill(found.begin(), found.end(), false);			
		}
		
		inline void update_quo_nodes(const int& i) {
			for (int j = 0; j < (int) temp.size(); j++) {
				in_set[temp[j]] = true;
			}
			
			int n = quo_nodes.size();
			for (int j = 0; j < n; j++) {
				while (j < n && in_set[quo_nodes[j]]) {
					swap(quo_nodes[j], quo_nodes[n-1]);
					quo_nodes.pop_back();
					n--;
				}
			}
			
			quo_nodes.push_back(i);
		}
		
		inline void all_neighbours(vector<int>& res, const int& i) {
			//a function similar to reach, should merge the two in the future
			res.clear();
			
			int w, k; queue<int> q;
			for (int j = 0; j < (int) graph[i].size(); j++) {
				w = graph[i][j];
				if (eliminated[w]) {
					q.push(w);
					found[w] = true;
				} else {
					res.push_back(w);
					found[w] = true;
				}
			}
			
			found[i] = true;
			while (!q.empty()) {
				k = q.front(); q.pop();

				for (int j = 0; j < (int) graph[k].size(); j++) {
					int w = graph[k][j];
					if (found[w]) continue;
					
					if (eliminated[w]) {
						q.push(w);
					} else {
						res.push_back(w);
						found[w] = true;
					}
				}
			}
			
			//cleaning up...
			//the line below can be optimized by keeping track of exactly what was put in.
			fill(found.begin(), found.end(), false);	
		}
		
		int union_size(int w, vector<int>& nbrs, int i) {
			all_neighbours(temp_nbrs, w);
			//cout << nbrs_w << endl;
						
			int res = 0;
			found[w] = found[i] = true;
			for (int j = 0; j < (int) temp_nbrs.size(); j++) {
				int v = temp_nbrs[j];
				if (!found[v] && !eliminated[v]) {
					res++;
					found[v] = true;
				}
			}
			
			for (int j = 0; j < (int) nbrs.size(); j++) {
				int v = nbrs[j];
				if (!found[v] && !eliminated[v]) {
					res++;
					found[v] = true;
				}
			}
			
			fill(found.begin(), found.end(), false);
			return res;
		}
		
		void scan() {
			for (int i = 0; i < (int) curr_nbrs.size(); i++) {
				int v = curr_nbrs[i];
				for (int j = 0; j < (int) graph[v].size(); j++) {
					int e = graph[v][j];
					if (degree[e] != -1) continue;
					if (L_approx[e] < 0) {
						all_neighbours(temp_nbrs, e);
						L_approx[e] = temp_nbrs.size();
					}
					L_approx[e] -= 1;
				}
			}
		}
		
		int approx_degree(int w) {
			//choice 1, number of nodes left
			int c1 = graph.size()-num_elim;
			
			//choice 2, degree of node previously + size of the clique formed on this step
			int Lw = curr_nbrs.size();
			int c2 = degree[w] + Lw;
			
			//choice 3, size of adjacency of ith node + size of clique formed + degrees of neighbouring quotient nodes.
			int Ai = 0;
			int Le = 0;
			for (int i = 0; i < (int) graph[w].size(); i++) {
				int v = graph[w][i];
				if (degree[v] != -1) {
					Ai++;
				} else {
					Le += L_approx[v];
				}
			}
			
			int c3 = Ai + Lw + Le;
			return min(c1, min(c2, c3));
		}
		
		void eliminate(const int& i) {
			// cout << " ============================================= " << endl;
			
			//1. temp <- common elements in Adj(i) and quo_nodes
			//2. reach <- reachable set of (v_i and quo_nodes)
			//3. v <- temp and v_i
			//4. quo_nodes <- {quo_nodes-temp} and {v}
			//5. Adj(v) <- reach
			//6. for every vertex w in reach:
			//	    Adj(w) <- {v} and Adj(w)-(temp and v_i)
			
			// cout << "on iter: " << i << endl;
			update_temp(i);      //1
			// cout << "done update" << endl;
			reachable(i);        //2
			// cout << "done reachable" << endl;
			update_quo_nodes(i);  //3 && 4
			// cout << "done quo" << endl;
			eliminated[i] = true;
			num_elim++;

			//5.
			//merges everything one level down. after merging, will try to
			//merge another level further if there is a successful merge
			int w, n = graph[i].size();
			for (int j = 0; j < n; j++) {
				w = graph[i][j];
				int m = graph[w].size();
				while (eliminated[w]) {
					if (m <= MERGE_CAP) {
						swap(graph[i][j], graph[i][n-1]);
						graph[i].pop_back();
						n--;

						// cout << graph[w] << endl;
						for (int k = 0; k < m; k++) {
							if (graph[w][k] != i) {
								graph[i].push_back(graph[w][k]);
								n++;
							}
						}
						// cout << "clearing " << w << endl;
						graph[w].clear();
						
						if (j >= n) break;
						w = graph[i][j];
						m = graph[w].size();
					} else {
						break;
					}
				}
			}
			
						
			//keep only uniques
			sort (graph[i].begin(), graph[i].end());
			vector<int>::iterator it = unique (graph[i].begin(), graph[i].end());
			graph[i].resize(distance(graph[i].begin(), it));
			
			all_neighbours(curr_nbrs, i);
			degree[i] = -1; //eliminated nodes are set to -1 for degree
			//end 5
			
			//6.
			//removing back references to current node in all previous quonodes.
			for (int j = 0; j < (int) temp.size(); j++) {
				in_set[temp[j]] = true;
				
				int m = graph[temp[j]].size();
				for (int k = 0; k < m; k++) {
					if (graph[temp[j]][k] == i) {
						swap(graph[temp[j]][k], graph[temp[j]][m-1]);
						graph[temp[j]].pop_back();
						break;
					}
				}
			}
			in_set[i] = true;
			
			
			//updating neighbours of i with new degrees
			//the three loops below must be done in order!
			for (int k = 0; k < (int) reach.size(); k++) {
				w = reach[k];
				if (eliminated[w]) {
					graph[w].clear();
					continue;
				}
				int m = graph[w].size();
				for (int j = 0; j < m; j++) {
					while (j < m && in_set[graph[w][j]]) {
						swap(graph[w][j], graph[w][m-1]);
						graph[w].pop_back();
						m--;
					}
				}
			}
			
			// cout << "updating degrees..." << endl;
			//update degrees. this part changes for md!
			scan();
			// cout << "done scan" << endl;
			for (int k = 0; k < (int) curr_nbrs.size(); k++) {
				w = curr_nbrs[k];
				degree[w] = approx_degree(w);
			}
			
			for (int k = 0; k < (int) reach.size(); k++) {
				w = reach[k];
				if (eliminated[w]) {
					graph[w].clear();
					continue;
				}
				graph[w].push_back(i);
			}
			
			for (int j = 0; j < (int) temp.size(); j++) {
				in_set[temp[j]] = false;
			}
			in_set[i] = false;
			
			//end 6
			
			// cout << "temp: " << temp << endl;
			// cout << "reach: " << reach << endl;
			// cout << "quo: " << quo_nodes << endl;
			// cout << "in_set: " << in_set << endl;
			// print();
		}
		
		void print() {
			for (int i = 0; i < (int) graph.size(); i++) {
				cout << i << ": " << graph[i] << '\n';
			}
			cout << "degree: " << degree << endl;
		}
		
		
};