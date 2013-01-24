
//overall todo list:
//	change variable names to make more sense.
//  documentation
//  remove explicit using namespace, etc.

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

class quotient_graph {
	typedef vector<vector<int> > adj_list;
	
	adj_list graph;
	vector<int> link, deg;
	vector<bool> in_set;
	vector<int> reach, quo_nodes, temp;
	
	const int MERGE_CAP = 3;
	public:
		void quotient_graph(adj_list g) : graph(g) {
			int n = g.size();
			link.resize(n, -1);
			deg.resize(n);
			in_set.resize(n, false);
			temp.reserve(n); 
			
			for (int i = 0; i < g.size(); i++) {
				deg[i] = g[i].size();
			}
		}
		
		inline void update_temp(const int& i) {
			temp.clear();
			for (int j = 0; j < graph[i].size(); j++) {
				in_set[graph[i][j]] = true;
			}
			
			for (int j = 0; j < quo_nodes.size(); j++) {
				if (in_set[quo_nodes[j]]) 
					temp.push_back(j);
			}
			
			for (int j = 0; j < graph[i].size(); j++) {
				in_set[graph[i][j]] = false;
			}
		}
		
		inline void reachable(const int& i) {
			reach.clear();
			for (int j = 0; j < quo_nodes.size(); j++) {
				in_set[quo_nodes[j]] = true;
			}
			
			for (int j = 0; j < graph[i].size(); j++) {
				if (!in_set[graph[i][j]]) {
					reach.push_back(graph[i][j]);
					in_set[graph[k][i]] = true;
				}
			}
			
			//update the links below!
			// int j = 0, k = i;
			// while (true) {
				// if (!in_set[graph[k][j]) {
					// reach.push_back(graph[k][j]);
					// in_set[graph[k][j]] = true;
				// }
				// j++;
				
				// this is pretty ugly code below.
				// refactor after doing other parts.
				// if (j >= graph[k].size() || graph[k][j] != -1) {
					// if (link[k] != -1) {
						// k = link[k];
						// j = 0;
					// }
				// } else {
					// break;
				// }
			// }
			
			//cleaning up...
			for (int j = 0; j < quo_nodes.size(); j++) {
				in_set[quo_nodes[j]] = false;
			}
			
			for (int j = 0; j < reach.size(); j++) {
				in_set[reach[j]] = false;
			}
			
		}
		
		inline void update_quo_node(const int& i) {
			for (int j = 0; j < temp.size(); j++) {
				in_set[temp[j]] = true;
			}
			
			for (int j = 0; j < quo_nodes.size(); j++) {
				while (j < quo_nodes.size() && in_set[quo_nodes[j]]) {
					swap(quo_nodes[j], quo_nodes[quo_nodes.size()-1]);
					quo_nodes.pop_back();
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
			//	    Adj(w) <- {v} and Adj(w)-(T and v_i)
			update_temp(i);      //1
			reachable(i);        //2
			update_quo_nodes();  //3 && 4
			
			//5.
			int w;
			for (int j = 0; j < graph[i].size(); j++) {
				w = graph[i][j];
				while (w < i) {
					swap(graph[i][j], graph[i][graph[i].size()-1]);
					graph[i].pop_back();
					if (graph[w].size() < CAP) {
						for (int k = 0; k < graph[w].size(); k++) {
							if (graph[w][k] != i) {
								graph[i].push_back(graph[w][k]);
							}
						}
						graph[i].clear();
					}
				}
			}
			//end 5
			
			//6.
			for (int j = 0; j < temp.size(); j++) {
				in_set[temp[j]] = true;
			}
			in_set[i] = true;
			
			for (int w = 0; w < reach.size(); w++) {
				for (int j = 0; j < graph[w].size(); j++) {
					while (j < graph[w].size() && in_set[graph[w][j]]) {
						swap(graph[w][j], graph[w][graph[w].size()-1]);
						graph[w].pop_back();
					}
				}
				graph[w].push_back(i);
			}
			//end 6
			
		}
		
		void print() {
			for (int i = 0; i < graph.size(); i++) {
				cout << i << ": " << graph[i] << endl;
			}
		}
		
		
};