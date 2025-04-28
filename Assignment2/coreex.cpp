#include <iostream>
#include <unordered_set>
#include <set>
#include <fstream>
using namespace std;


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
void BronKerboschPivot(vector<int> &P, vector<int> &R, vector<int> &X, vector<set<int>> &adj, int n, int h, vector<vector<int>>& hCliques, vector<int>& cliqueDegrees, vector<vector<int>>& vCliques) {
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
            vCliques[x].push_back(hCliques.size() - 1);
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
        BronKerboschPivot(intersectP, R, intersectX, adj, n, h, hCliques, cliqueDegrees, vCliques);
        R.pop_back();
        X.push_back(i);
        it = P.erase(it);
    }
}


struct Edge {
    int dest, revIdx;
    double capacity;
};

void addEdge(vector<vector<Edge>>& graph, int src, int dest, double capacity) {
    graph[src].push_back({dest, (int)graph[dest].size(), capacity});
    graph[dest].push_back({src, (int)graph[src].size() - 1, 0.0});
}

double eps = 1e-9;
double INF = 1e14;


double coreDecomposition(int n, int h, vector<set<int>> &adj, vector<int>& cliqueDegrees, vector<vector<int>>& hCliques, vector<vector<int>> &vCliques, vector<int>& kCoreVertices) {
    int maxDegree = 0;
    for (int i = 0; i < n; i++) {
        maxDegree = max(maxDegree, cliqueDegrees[i]);
    }
    vector<int> bin(maxDegree + 1, 0);
    for (int i = 0; i < n; i++) {
        bin[cliqueDegrees[i]]++;
    }
    int start = 0;
    for (int i = 0; i <= maxDegree; i++) {
        int num = bin[i];
        bin[i] = start;
        start += num;
    }
    vector<int> pos(n);
    vector<int> vert(n);
    for (int i = 0; i < n; i++) {
        pos[i] = bin[cliqueDegrees[i]];
        vert[pos[i]] = i;
        bin[cliqueDegrees[i]]++;
    }
    for (int i = maxDegree; i >= 1; i--) {
        bin[i] = bin[i - 1];
    }
    bin[0] = 1;
    int completed = 0;
    double hcDensity = ((double) hCliques.size()) / n;
    for (int i = 0; i < n; i++) kCoreVertices.push_back(i);
    vector<int> vertexRemoved(n, 0);

    vector<int> R;
    vector<vector<int>> hCliquesBK;
    vector<int> cliqueDegreesBK(n);
    vector<vector<int>> vCliquesTemp(n);
    int num_cliques = hCliques.size();
    vector<int> cliqueRemoved(hCliques.size(), 0);
    for (int i = 0; i < n; i++) {
        int minVal = INT_MAX;
        int minIdx = 0;
        for (int j = 0; j < n; j++) {
            if (vertexRemoved[j] == 0) {
                if (cliqueDegrees[j] < minVal) {
                    minVal = cliqueDegrees[j];
                    minIdx = j;
                }
            }
        }
        int v = minIdx;
        for (auto cliqueIdx: vCliques[v]) {
            if (cliqueRemoved[cliqueIdx] == 0) {
                for (auto u: hCliques[cliqueIdx]) {
                    if (cliqueDegrees[u] > cliqueDegrees[v]) cliqueDegrees[u]--;
                }
            }
        }
        completed++;
        if (vertexRemoved[v] == 1) cout << "removing duplicate vertex " << v << endl;
        vertexRemoved[v] = 1;
        for (int clique = 0; clique < (int)vCliques[v].size(); clique++) {
            if (cliqueRemoved[vCliques[v][clique]] == 0) {
                cliqueRemoved[vCliques[v][clique]] = 1;
                num_cliques--;
            }
        }
        double val =  ((double)num_cliques / (n - completed));
        if (val > hcDensity) {
            kCoreVertices.clear();
            for (int m = 0; m < n; m++) {
                if (vertexRemoved[m] == 0) kCoreVertices.push_back(m);
            }
            hcDensity = val;
        }
        
    }
    return hcDensity;
}


void dfs(int i, vector<set<int>>& adj, vector<bool>& visited, vector<int>& connectedComponent) {
    visited[i] = true;
    connectedComponent.push_back(i);
    for (auto x: adj[i]) {
        if (!visited[x]) dfs(x, adj, visited, connectedComponent);
    }
}

double flowDfs(vector<vector<Edge>>& graph, int i, int t, double flow, vector<bool>& visited) {
    if (i == t) return flow;
    visited[i] = true;
    for (Edge& e : graph[i]) {
        if (!visited[e.dest] && e.capacity > eps) {
            double bottleneck = flowDfs(graph, e.dest, t, min(flow, e.capacity), visited);
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
        double bottleneck = flowDfs(graph, s, t, INF, visited);
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



double connectedComponents(int n, vector<set<int>> &adj, vector<int>& inKCore, vector<vector<int>>& hCliques, vector<vector<int>>& vCliques, vector<vector<int>>& allConnectedComponents, int &numVertices) {
    vector<bool> visited(n, false);
    vector<int> connectedComponent;
    double kdoubleprime = 0.0;
    for (int i = 0; i < n; i++) {
        if (inKCore[i]) {
            if (!visited[i]) {
                connectedComponent.clear();
                int count_cliques = 0;
                dfs(i, adj, visited, connectedComponent);
                vector<int> cliqueChecked(hCliques.size(), 0);
                allConnectedComponents.push_back(connectedComponent);
                for (auto x: connectedComponent) {
                    for (auto clique: vCliques[x]) {
                        int flag = 0;
                        if (cliqueChecked[clique] == 0) {
                            for (auto u: hCliques[clique]) {
                                if (inKCore[u] == 0) {
                                    flag = 1;
                                    break;
                                }
                            }
                            if (flag == 0) count_cliques++;
                            cliqueChecked[clique] = 1;
                        }
                    }
                }
                if (((double) count_cliques / (double) connectedComponent.size()) > kdoubleprime) {
                    kdoubleprime = ((double) count_cliques / (double) connectedComponent.size());
                    numVertices = connectedComponent.size();
                }
            }
        }
    }
    return kdoubleprime;
}

void getConnectedComponents(int n, vector<set<int>>& adj, vector<int>& inKCore, vector<vector<int>>& allConnectedComponents) {
    vector<bool> visited(n, false);
    for (int i = 0; i < n; i++) {
        if (inKCore[i]) {
            if (!visited[i]) {
                vector<int> connectedComponent;
                dfs(i, adj, visited, connectedComponent);
                allConnectedComponents.push_back(connectedComponent);
            }
        }
    }
}


void buildFlowNetwork(int n, vector<set<int>>& adj, vector<vector<Edge>>& graph, vector<vector<int>>& deltaCC, vector<int>& cliqueDegrees, double alpha, int h) {
    int numCliques = deltaCC.size();
    int s = n;
    int t = n + 1;
    int totalVertices = n + 2 + numCliques;
    graph.resize(totalVertices);
    for (int i = 0; i < n; i++) {
        addEdge(graph, s, i, (double) cliqueDegrees[i]);
        addEdge(graph, i, t, alpha * (double) h);
       
    }
    for (int i = 0; i < numCliques; i++) {
        for (auto x: deltaCC[i]) {
            addEdge(graph, (n + 2 + i), x, (double) INF);
        }
    }
    for (int i = 0; i < numCliques; i++) {
        for (int v = 0; v < n; v++) {
            int flag = 1;
            for (auto x: deltaCC[i]) {
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
}


void coreExact(int n, int h, vector<set<int>>& adj, vector<vector<int>>& hCliques, vector<int>& cliqueDegrees, vector<vector<int>>& vCliques, vector<vector<int>>& delta, vector<vector<int>>& deltaMap, vector<int>& D) {
    vector<int> tempDegrees;
    for (auto x: cliqueDegrees) tempDegrees.push_back(x);
    vector<int> kCoreVertices;
    double hcDensity = coreDecomposition(n, h, adj, tempDegrees, hCliques, vCliques, kCoreVertices);
    int kprime = ceil(hcDensity);
    vector<set<int>> kCoreAdj(n);
    vector<int> inKCore(n, 0);
    int coreVertices = kCoreVertices.size();
    for (int i = 0; i < n; i++) {
        if (tempDegrees[i] >= kprime) {
            inKCore[i] = 1;
            for (auto x: adj[i]) {
                if (tempDegrees[x] >= kprime) {
                    kCoreAdj[i].insert(x);
                    
                } 
            }
        }
    }
    vector<vector<int>> allConnectedComponents;
    int numVertices = 0;
    double ccDensity = connectedComponents(n, kCoreAdj, inKCore, hCliques, vCliques, allConnectedComponents, numVertices);
    int kdoubleprime;
    if (ceil(ccDensity) > kprime) {
        kdoubleprime = ceil(ccDensity);
        vector<set<int>> kdpCoreAdj(n);
        vector<int> inKdpCore(n, 0);
        for (int i = 0; i < n; i++) {
            if (tempDegrees[i] >= kdoubleprime) {
                inKdpCore[i] = 1;
                for (auto x: adj[i]) {
                    if (tempDegrees[x] >= kdoubleprime) kdpCoreAdj[i].insert(x);
                }
            }
        }
        coreVertices = numVertices;
        allConnectedComponents.clear();
        getConnectedComponents(n, kdpCoreAdj, inKdpCore, allConnectedComponents);
    } else kdoubleprime = kprime;
    if (allConnectedComponents.size() == 0) allConnectedComponents.push_back(kCoreVertices);
    int kmax = 0;
    for (auto x: tempDegrees) kmax = max(kmax, x);
    double l = kdoubleprime;
    double u = (double)coreVertices + 1.0;
    cout << "Inital Parameters:\n";
    cout << "l: " << l << " u: " << u << endl;
    if (allConnectedComponents.size() > 0) {
        for (auto x: allConnectedComponents[0]) {
            D.push_back(x);
        }
    }
    vector<int> U; 
    for (int i = 0; i < allConnectedComponents.size(); i++) {
        if (l > kdoubleprime) {
            vector<int> newC;
            for (auto x: allConnectedComponents[i]) {
                if (tempDegrees[x] >= ceil(l)) newC.push_back(x);
            }
            allConnectedComponents[i] = newC;
        }
        vector<int> inCC(n, 0);
        for (auto x: allConnectedComponents[i]) inCC[x] = 1;

        vector<int> ccVertices; 
        vector<int> oldToNew(n, -1); 
        vector<int> newToOld;      

        for (int v = 0; v < n; v++) {
            if (inCC[v] == 1) {
                oldToNew[v] = ccVertices.size();
                newToOld.push_back(v);
                ccVertices.push_back(v);
            }
        }
        int ccSize = ccVertices.size();

        vector<set<int>> ccAdj(ccSize);
        for (int v : ccVertices) {
            int newV = oldToNew[v];
            for (auto x : adj[v]) {
                if (inCC[x] == 1) {
                    int newX = oldToNew[x];
                    ccAdj[newV].insert(newX);
                }
            }
        }

        vector<vector<int>> ccDelta;
        vector<int> deltaChecked(delta.size(), 0);
        for (auto x : allConnectedComponents[i]) {
            for (auto clique : deltaMap[x]) {
                if (deltaChecked[clique] == 0) {
                    int flag = 0;
                    for (auto u : delta[clique]) {
                        if (inCC[u] == 0) {
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0) {
                        vector<int> temp;
                        for (auto x: delta[clique]) temp.push_back(oldToNew[x]);
                        ccDelta.push_back(temp);
                    }
                    deltaChecked[clique] = 1;
                }
            }
        }
        vector<vector<Edge>> graph;
        vector<int> cliqueChecked(hCliques.size(), 0);
        vector<int> newDegs(ccAdj.size(), 0);
        for (auto x: allConnectedComponents[i]) {
            for (auto clique: vCliques[x]) {
                if (cliqueChecked[clique] == 0) {
                    int flag = 0;
                    for (auto v: hCliques[clique]) {
                        if (inCC[v] == 0) {
                            flag = 1;
                            break;
                        }
                    }
                    if (flag == 0) {
                        for (auto v: ccVertices) newDegs[oldToNew[v]] = cliqueDegrees[v];
                    }

                    cliqueChecked[clique] = 1;
                }
            }
        }
        buildFlowNetwork(ccAdj.size(), ccAdj, graph, ccDelta, newDegs, (double) l, h);
        vector<int> S;
        minCut(graph, S, ccAdj.size());
        if (S.empty()) continue;
        vector<int> inU(n, 0);
        for (auto x: S) {
            if (x < ccAdj.size()) {
                U.push_back(newToOld[x]);
                inU[newToOld[x]] = 1;
            }
        } 
        if (!D.empty()) {
            int countUCliques = 0;
            for (auto x: U) {
                for (auto clique: vCliques[x]) {
                    if (cliqueChecked[clique] == 0) {
                        int flag = 0;
                        for (auto y: hCliques[clique]) {
                            if (inU[y] == 0) {
                                flag = 1;
                                break;
                            }
                        }
                        if (flag == 0) countUCliques++;
                        cliqueChecked[clique] = 1;
                    }
                }
            }
            fill(cliqueChecked.begin(), cliqueChecked.end(), 0);
            int countDCliques = 0;
            for (auto x: D) {
                for (auto clique: vCliques[x]) {
                    if (cliqueChecked[clique] == 0) {
                        int flag = 0;
                        for (auto y: hCliques[clique]) {
                            if (inU[y] == 0) {
                                flag = 1;
                                break;
                            }
                        }
                        if (flag == 0) countDCliques++;
                        cliqueChecked[clique] = 1;
                    }
                }
            }
            double uDensity = ((double) countUCliques / (double) U.size());
            double dDensity = ((double) countDCliques / (double) D.size());
        
            if (uDensity > dDensity) {
                D.clear();
                for (auto x: U) D.push_back(x);
            }
        } else {
            for (auto x: U) D.push_back(x);
        }
        cout << "l is: " << l << " and u is: " << u << endl;
        cout << "val is: " << (1.0/((double)allConnectedComponents[i].size() * (allConnectedComponents[i].size() - 1))) << endl;
        while ((u - l) >= (1.0/((double)allConnectedComponents[i].size() * (allConnectedComponents[i].size() - 1)))) {
            double alpha = (l + u) / 2.0;
            cout << "l is: " << l << " and u is: " << u << endl;
            cout << "alpha is: " << alpha << endl;
            cout << "val is: " << (1.0/((double)allConnectedComponents[i].size() * (allConnectedComponents[i].size() - 1))) << endl;
            graph.clear();
            buildFlowNetwork(ccAdj.size(), ccAdj, graph, ccDelta, newDegs, (double) alpha, h);
            S.clear();
            minCut(graph, S, ccAdj.size());
            if (S.size() == 1 && S[0] == ccAdj.size()) {
                u = alpha;
            } else {
                if (alpha > ceil(l)) {
                    fill(inCC.begin(), inCC.end(), 0);
                    vector<int> newC;
                    for (auto x: allConnectedComponents[i]) {
                        if (cliqueDegrees[x] >= ceil(l)) {
                            newC.push_back(x);
                            inCC[x] = 1;
                        }
                        ccAdj.clear();
                    }
                    allConnectedComponents[i] = newC;
                    ccVertices.clear();
                    fill(oldToNew.begin(), oldToNew.end(), -1);
                    newToOld.clear();
                    
                    for (int v = 0; v < n; v++) {
                        if (inCC[v] == 1) {
                            oldToNew[v] = ccVertices.size();
                            newToOld.push_back(v);
                            ccVertices.push_back(v);
                        }
                    }
                    ccSize = ccVertices.size();
                    ccAdj.clear();
                    ccAdj.resize(ccSize);
                    for (int v : ccVertices) {
                        int newV = oldToNew[v];
                        for (auto x : adj[v]) {
                            if (inCC[x] == 1) {
                                int newX = oldToNew[x];
                                ccAdj[newV].insert(newX);
                            }
                        }
                    }
                    
                    ccDelta.clear();
                    fill(deltaChecked.begin(), deltaChecked.end(), 0);
                    for (auto x : allConnectedComponents[i]) {
                        for (auto clique : deltaMap[x]) {
                            if (deltaChecked[clique] == 0) {
                                int flag = 0;
                                for (auto u : delta[clique]) {
                                    if (inCC[u] == 0) {
                                        flag = 1;
                                        break;
                                    }
                                }
                                if (flag == 0) {
                                    vector<int> temp;
                                    for (auto x: delta[clique]) temp.push_back(oldToNew[x]);
                                    ccDelta.push_back(temp);
                                }
                                deltaChecked[clique] = 1;
                            }
                        }
                    }

                    fill(cliqueChecked.begin(), cliqueChecked.end(), 0);
                    fill(newDegs.begin(), newDegs.end(), 0);
                    for (auto x: allConnectedComponents[i]) {
                        for (auto clique: vCliques[x]) {
                            if (cliqueChecked[clique] == 0) {
                                int flag = 0;
                                for (auto v: hCliques[clique]) {
                                    if (inCC[v] == 0) {
                                        flag = 1;
                                        break;
                                    }
                                }
                                if (flag == 0) {
                                    for (auto v: ccVertices) newDegs[oldToNew[v]] = cliqueDegrees[v];
                                }

                                cliqueChecked[clique] = 1;
                            }
                        }
                    }
                }
                l = alpha;
                U.clear();
                for (auto x: S) {
                    if (x < ccAdj.size()) {
                        inU[newToOld[x]] = 1;
                        U.push_back(newToOld[x]);
                    }

                }
            }
        }
        if (!D.empty()) {
            int countUCliques = 0;
            for (auto x: U) {
                for (auto clique: vCliques[x]) {
                    if (cliqueChecked[clique] == 0) {
                        int flag = 0;
                        for (auto y: hCliques[clique]) {
                            if (inU[y] == 0) {
                                flag = 1;
                                break;
                            }
                        }
                        if (flag == 0) countUCliques++;
                        cliqueChecked[clique] = 1;
                    }
                }
            }
            fill(cliqueChecked.begin(), cliqueChecked.end(), 0);
            int countDCliques = 0;
            for (auto x: D) {
                for (auto clique: vCliques[x]) {
                    if (cliqueChecked[clique] == 0) {
                        int flag = 0;
                        for (auto y: hCliques[clique]) {
                            if (inU[y] == 0) {
                                flag = 1;
                                break;
                            }
                        }
                        if (flag == 0) countDCliques++;
                        cliqueChecked[clique] = 1;
                    }
                }
            }
            double uDensity = ((double) countUCliques / (double) U.size());
            double dDensity = ((double) countDCliques / (double) D.size());
        
            if (uDensity > dDensity) {
                D.clear();
                for (auto x: U) D.push_back(x);
            }
        } else {
            for (auto x: U) D.push_back(x);
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


    vector<int> degeneracyOrder(n, -1);
    int degeneracy = getDegeneracy(adj, degeneracyOrder, n);
    vector<int> R;
    unordered_set<int> completed;
    vector<vector<int>> hCliques;
    vector<int> cliqueDegrees(n);
    vector<vector<int>> vCliques(n);
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
        BronKerboschPivot(P, R, X, adj, n, h, hCliques, cliqueDegrees, vCliques);
        R.pop_back();
        completed.insert(i);
    }
    cout << "Found " << hCliques.size() << " hCliques\n";
    R.clear();
    completed.clear();
    vector<vector<int>> delta;
    vector<int> temp(n);
    vector<vector<int>> deltaMap(n);
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
        BronKerboschPivot(P, R, X, adj, n, h - 1, delta, temp, deltaMap);
        R.pop_back();
        completed.insert(i);
    }
    cout << "Found " << delta.size() << " h - 1 Cliques\n";
    vector<int> D;
    coreExact(n, h, adj, hCliques, cliqueDegrees, vCliques, delta, deltaMap, D);
    cout << D.size() << endl;
    vector<int> inD(n, 0);
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
