#include <iostream>
#include <unordered_set>
#include <set>
#include <fstream>
#include <vector>   
#include <set>     
#include <algorithm>
using namespace std;

struct Edge {
    int dest, revIdx;
    double capacity;
};

double eps = 1e-14;
double INF = 1e14;


int getDegeneracy(vector<set<int>> &adj, vector<int>& degeneracyOrder, int n) {
    vector<int> degrees(n);
    for (int i = 0; i < n; i++) {
        degrees[i] = adj[i].size();
    }
    vector<vector<int>> D(n);
    for (int i = 0; i < n; i++) {
        D[degrees[i]].push_back(i);
    }
    int i = 0; 
    int num = 0;
    int degeneracy = 0;
    vector<bool> removed(n, false);
    while (num < n) {
        while (D[i].empty()) i++;
        int x = D[i].front();
        D[i].erase(D[i].begin());
        degeneracyOrder[num++] = x;
        degeneracy = max(degeneracy, i);
        removed[x] = true;
        for (auto y: adj[x]) {
            if (!removed[y]) {
                int oldDegree = degrees[y];
                degrees[y]--;
                auto it = find(D[oldDegree].begin(), D[oldDegree].end(), y);
                if (it!= D[oldDegree].end()) {
                    D[oldDegree].erase(it);
                }
                
                D[degrees[y]].push_back(y);
            }
        }
        i = 0;
    }
    return degeneracy;
}
void BronKerboschPivot(vector<int> &P, vector<int> &R, vector<int> &X, vector<set<int>> &adj, int n, int h, vector<vector<int>>& hCliques, vector<int>& cliqueDegrees) {
    if ((int)R.size() == h) {
        vector<int> inR(n, 0);
        int flag = 0;
        for (auto x: R) {
            if (inR[x] == 1) {
                flag = 1;
                break;
            }
            inR[x] = 1;
        }
        if (flag == 1) return;
        hCliques.push_back(R);
        for (auto x: R) {
            cliqueDegrees[x]++;
        }
        return;
    }
    for (auto it = P.begin(); it != P.end();) {
        int i = *it;
        vector<int> intersectP;
        vector<int> intersectX;
        for (auto y: P) {
            if (adj[i].find(y) != adj[i].end()) intersectP.push_back(y);
        }
        for (auto y: X) {
            if (adj[i].find(y) != adj[i].end()) intersectX.push_back(y);
        }

        R.push_back(i);
        BronKerboschPivot(intersectP, R, intersectX, adj, n, h, hCliques, cliqueDegrees);
        R.pop_back();
        X.push_back(i);
        it = P.erase(it);
    }
}


void addEdge(vector<vector<Edge>>& graph, int src, int dest, double capacity) {
    graph[src].push_back({dest, (int)graph[dest].size(), capacity});
    graph[dest].push_back({src, (int)graph[src].size() - 1, 0.0});
}

double dfs(vector<vector<Edge>>& graph, int i, int t, double flow, vector<bool>& visited) {
    if (i == t) return flow;
    visited[i] = true;
    for (Edge& e : graph[i]) {
        if (!visited[e.dest] && e.capacity > eps) {
            double bottleneck = dfs(graph, e.dest, t, min(flow, e.capacity), visited);
            if (bottleneck > eps) {
                e.capacity -= bottleneck;
                graph[e.dest][e.revIdx].capacity += bottleneck;
                return bottleneck;
            }
        }
    }
    return 0.0;
}

void maxFlow(vector<vector<Edge>>& graph, int s, int t) {
    double flow = 0.0;
    int n = graph.size();
    while (true) {
        vector<bool> visited(n, false);
        double bottleneck = dfs(graph, s, t, INF, visited);
        if (bottleneck < eps) break;
        flow += bottleneck;
    }
    cout << "Max Flow: " << flow << endl;
}

void minCut(vector<vector<Edge>>& graph, vector<int>& S, int n) {
    maxFlow(graph, n, n + 1);
    vector<bool> visited(graph.size(), false);
    stack<int> st;
    st.push(n);
    visited[n] = true;
    while (!st.empty()) {
        int v = st.top();
        st.pop();
        S.push_back(v);
        for (Edge& e: graph[v]) {
            if (!visited[e.dest] && e.capacity > eps) {
                visited[e.dest] = true;
                st.push(e.dest);
            }
        }
    }
}

void exact(int n, vector<set<int> >& adj, int h, vector<vector<int>>& hCliques, vector<int>& cliqueDegrees, vector<vector<int>>& delta, vector<int>& D, unordered_map<int, int>& mapping, unordered_map<int, int>& revMap) {
    double l = 0.0;
    double u = 0.0;
    for (int i = 0; i < n; i++) {
        u = max(u, (double) cliqueDegrees[i]);
    }
    cout << "inital parameters: " << endl;
    cout << "l = " << l << " u = " << u << endl;
    int numCliques = delta.size();
    double val = (1.0/((double)n * (n - 1)));

    cout << "val is: " << val << endl;
    while ((u - l) >= val) {
        cout << "l is: " << l << " u is: " << u << endl;
        double alpha = (l + u) / 2.0;
        cout << "Alpha is: " << alpha << endl;
        int s = n;
        int t = n + 1;
        int totalVertices = n + 2 + numCliques;
        vector<vector<Edge>> graph(totalVertices);
        for (int i = 0; i < n; i++) {
            addEdge(graph, s, i, (double) cliqueDegrees[i]);
            addEdge(graph, i, t, alpha * (double) h);
        }
        for (int i = 0; i < numCliques; i++) {
            for (auto x: delta[i]) {
                addEdge(graph, (n + 2 + i), x, (double) INF);
            }
        }
        for (int i = 0; i < numCliques; i++) {
            for (int v = 0; v < n; v++) {
                int flag = 1;
                for (auto x: delta[i]) {
                    if (adj[x].find(v) == adj[x].end()) {
                        flag = 0;
                        break;
                    }
                }
                if (flag == 1) {
                    addEdge(graph, v, (n + 2 + i), 1.0);
                }
            }
        }
        vector<int> S;
        minCut(graph, S, n);
        cout << endl;
        if (S.size() == 1 && S[0] == n) {
            u = alpha;
        } else {
            l = alpha;
            D.clear();
            for (auto x: S) {
                if (x < n) D.push_back(x);
            }
        }
    }


}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }
    auto start = chrono::high_resolution_clock::now();
    string filename = argv[1];
    vector<set<int>> adj;
    int n, m, h;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    infile >> h;       
    infile >> n >> m; 
    unordered_map<int, int> map_input_to_internal;
    unordered_map<int, int> rev_map;
    int next_id = 0;

    adj.resize(n);  

    int u, v;
    for (int i = 0; i < m; ++i) {
        infile >> u >> v;
        if (u != v) {
            if (map_input_to_internal.find(u) == map_input_to_internal.end()) {
                map_input_to_internal[u] = next_id++;
            }
            if (map_input_to_internal.find(v) == map_input_to_internal.end()) {
                map_input_to_internal[v] = next_id++;
            }

            int u_internal = map_input_to_internal[u];
            int v_internal = map_input_to_internal[v];

            adj[u_internal].insert(v_internal);
            adj[v_internal].insert(u_internal);
        }
    }
    cout << "n is: " << n << endl;
    cout << "m is: " << m << endl;

    infile.close();

    vector<int> degeneracyOrder(n, -1);
    int degeneracy = getDegeneracy(adj, degeneracyOrder, n);
    vector<int> R;
    unordered_set<int> completed;
    vector<vector<int>> hCliques;
    vector<int> cliqueDegrees(n);
    for (int t = 0; t < n; t++) {
        int i = degeneracyOrder[t];
        vector<int> P;
        vector<int> X;
        for (auto x: adj[i]) {
            if (completed.find(x) != completed.end()) {
                X.push_back(x);
            } else P.push_back(x);
        }
        R.push_back(i);
        BronKerboschPivot(P, R, X, adj, n, h, hCliques, cliqueDegrees);
        R.pop_back();
        completed.insert(i);
    }
    cout << "Found " << hCliques.size() << " hCliques:\n";
    R.clear();
    completed.clear();
    vector<vector<int>> delta;
    vector<int> temp(n);
    for (int t = 0; t < n; t++) {
        int i = degeneracyOrder[t];
        vector<int> P;
        vector<int> X;
        for (auto x: adj[i]) {
            if (completed.find(x) != completed.end()) {
                X.push_back(x);
            } else P.push_back(x);
        }
        R.push_back(i);
        BronKerboschPivot(P, R, X, adj, n, h - 1, delta, temp);
        R.pop_back();
        completed.insert(i);
    }
    cout << "Found " << delta.size() << " h - 1 Cliques\n";
    vector<int> D;
    unordered_map<int, int> mapping;
    unordered_map<int, int> revMap;
    exact(n, adj, h, hCliques, cliqueDegrees, delta, D, mapping, revMap);
    cout << D.size() << endl;
    vector<int> inD(n, 0);
    cout << "D: " << endl;
    for (auto x: D) {
        inD[x] = 1;
        cout << x << " ";
    }
    cout << endl;
    int count_cliques = 0;
    for (int i = 0; i < hCliques.size(); i++) {
        int flag = 0;
        for (auto x: hCliques[i]) {
            if (inD[x] == 0) {
                flag = 1;
                break;
            }
        }
        if (flag == 0) count_cliques++;
    }
    auto end = chrono::high_resolution_clock::now();
    cout << "CDS Density: " << (((double) count_cliques) / (double) D.size()) << endl;
    cout << "Number of cliques: " << count_cliques << endl;
    cout << "Number of vertices: " << D.size() << endl;
    chrono::duration<double> duration = end - start;
    cout << "Time taken by algorithm: " << duration.count() << " seconds" << endl;
    return 0;
}