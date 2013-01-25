
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

const int MERGE_CAP = 3;
class quotient_graph {
	typedef vector<vector<int> > adj_list;
	
	adj_list graph;
	vector<bool> in_set;
	vector<int> reach, quo_nodes, temp;
	
	public:
		quotient_graph(adj_list& g) {
			int n = g.size();
			in_set.resize(n, false);
			temp.reserve(n); 
			reach.reserve(n);
			quo_nodes.reserve(n);
			graph.resize(n);
			
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < (int) g[i].size(); j++) {
					int v = g[i][j];
					if (v == i) continue;
					graph[i].push_back(v);
					graph[v].push_back(i);
				}
			}
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
				if (graph[i][j] < i) {
					q.push(graph[i][j]);
				}
				in_set[graph[i][j]] = true;
			}
		
			while (!q.empty()) {
				k = q.front(); q.pop();

				for (int j = 0; j < (int) graph[k].size(); j++) {
					w = graph[k][j];
					if (w == i) {
						swap(graph[k][j], graph[k][graph[k].size()-1]);
						graph[k].pop_back();
						reach.push_back(graph[k][j]);
						continue;
					}
					
					while (j < (int) graph[k].size() && in_set[graph[k][j]]) {
						swap(graph[k][j], graph[k][graph[k].size()-1]);
						graph[k].pop_back();
					}
					
					if (w < k) {
						q.push(w);
					} else {
						reach.push_back(w);
						in_set[w] = true;
					}
				}
			}			
			
			//cleaning up...
			for (int j = 0; j < (int) quo_nodes.size(); j++) {
				in_set[quo_nodes[j]] = false;
			}
			
			for (int j = 0; j < (int) reach.size(); j++) {
				in_set[reach[j]] = false;
			}
			
			for (int j = 0; j < (int) graph[i].size(); j++) {
				in_set[graph[i][j]] = false;
			}
			
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
		
		void eliminate(const int& i) {
			print();
			cout << " ============================================= " << endl;
			
			//1. temp <- common elements in Adj(i) and quo_nodes
			//2. reach <- reachable set of (v_i and quo_nodes)
			//3. v <- temp and v_i
			//4. quo_nodes <- {quo_nodes-temp} and {v}
			//5. Adj(v) <- reach
			//6. for every vertex w in reach:
			//	    Adj(w) <- {v} and Adj(w)-(temp and v_i)
			update_temp(i);      //1
			cout << "done update" << endl;
			reachable(i);        //2
			cout << "done reachable" << endl;
			update_quo_nodes(i);  //3 && 4
			cout << "done quo" << endl;
			
			//5.
			int w, n = graph[i].size();
			for (int j = 0; j < n; j++) {
				w = graph[i][j];
				int m = graph[w].size();
				while (w < i) {
					if (m <= MERGE_CAP) {
						swap(graph[i][j], graph[i][n-1]);
						graph[i].pop_back();
						n--;
						
						for (int k = 0; k < m; k++) {
							if (graph[w][k] != i) {
								graph[i].push_back(graph[w][k]);
								n++;
							}
						}
						graph[w].clear();
						
						w = graph[i][j];
					} else {
						break;
					}
				}
			}
			
			//keep only uniques
			sort (graph[i].begin(), graph[i].end());
			vector<int>::iterator it = unique (graph[i].begin(), graph[i].end());
			graph[i].resize(distance(graph[i].begin(), it));
			//end 5
			
			//6.
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
			
			for (int k = 0; k < (int) reach.size(); k++) {
				w = reach[k];
				int m = graph[w].size();
				for (int j = 0; j < m; j++) {
					cout << graph[w][j] << " ";
					while (j < m && in_set[graph[w][j]]) {
						swap(graph[w][j], graph[w][m-1]);
						graph[w].pop_back();
						m--;
					}
				}
				cout << endl;
				graph[w].push_back(i);
			}
			
			for (int j = 0; j < (int) temp.size(); j++) {
				in_set[temp[j]] = false;
			}
			in_set[i] = false;
			
			//end 6
			
			cout << "temp: " << temp << endl;
			cout << "reach: " << reach << endl;
			cout << "quo: " << quo_nodes << endl;
			cout << "in_set: " << in_set << endl;
		}
		
		void print() {
			for (int i = 0; i < (int) graph.size(); i++) {
				cout << i << ": " << graph[i] << endl;
			}
		}
		
		
};