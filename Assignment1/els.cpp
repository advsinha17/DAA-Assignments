#include <iostream>
#include <vector>
#include <list>
#include <set>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <unordered_set>
using namespace std;

long long count_cliques = 0;
long long max_size = 0;


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
void BronKerboschPivot(vector<int> &P, vector<int> &R, vector<int> &X, vector<set<int>> &adj, int n, unordered_map<int, int> &distribution) {
    if (P.empty() && X.empty()) {
        distribution[(int)R.size()]++;
        count_cliques++;
        max_size = max(max_size, (long long)R.size());
        return;
    }
    int pivot = -1;
    int pivotSize = -1;
    for (auto i: P) {
        int count = 0;
        for (auto y: P) {
            if (adj[i].find(y) != adj[i].end()) count++;
        }
        if (count > pivotSize) {
            pivot = i;
            pivotSize = count;
        }
        
    }
    for (auto i: X) {
        int count = 0;
        for (auto y: P) {
            if (adj[i].find(y) != adj[i].end()) count++;
        }
        if (count > pivotSize) {
            pivot = i;
            pivotSize = count;
        }
    }
    for (auto it = P.begin(); it != P.end();) {
        int i = *it;  

        if (adj[pivot].find(i) != adj[pivot].end()) {
            it++;
            continue;
        }
        vector<int> intersectP;
        vector<int> intersectX;
        for (auto y: P) {
            if (adj[i].find(y) != adj[i].end()) intersectP.push_back(y);
        }
        for (auto y: X) {
            if (adj[i].find(y) != adj[i].end()) intersectX.push_back(y);
        }

        R.push_back(i);
        BronKerboschPivot(intersectP, R, intersectX, adj, n, distribution);
        R.pop_back();
        X.push_back(i);
        it = P.erase(it);
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Correct Usage: %s <dataFilePath>\n", argv[0]);
        exit(1);
    }
    auto start = chrono::high_resolution_clock::now();

    string filename = argv[1];
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening file" << endl;
        return 1;
    }

    string line;
    int n = 0, m = 0;
    vector<set<int>> adj;
    vector<int> degrees;
    unordered_map<int, int> vertexMap; 
    int nextIndex = 0;       

    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            if (line.find("Nodes:") != string::npos) {
                size_t pos = line.find("Nodes:");
                if (pos != string::npos) {
                    n = stoi(line.substr(pos + 7));
                }
            }
            if (line.find("Edges:") != string::npos) {
                size_t pos = line.find("Edges:");
                if (pos != string::npos) {
                    m = stoi(line.substr(pos + 7));
                }
            }
            continue;
        }

        if (n == 0 || m == 0) {
            continue;
        }

        istringstream iss(line);
        int u, v;
        if (!(iss >> u >> v)) {
            continue;
        }

        if (vertexMap.find(u) == vertexMap.end()) {
            vertexMap[u] = nextIndex;
            nextIndex++;
        }
        if (vertexMap.find(v) == vertexMap.end()) {
            vertexMap[v] = nextIndex;
            nextIndex++;
        }

        if (adj.size() < n) {
            adj.resize(n);
            degrees.resize(n, 0);
        }

        int mappedU = vertexMap[u];
        int mappedV = vertexMap[v];

        adj[mappedU].insert(mappedV);
        adj[mappedV].insert(mappedU);
    }

    infile.close();
    vector<int> degeneracyOrder(n, -1);
    int degeneracy = getDegeneracy(adj, degeneracyOrder, n);
    unordered_map<int, int> distribution;


    // BronKerboschDegeneracy Algorithm
    vector<int> R;
    unordered_set<int> completed;
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
        BronKerboschPivot(P, R, X, adj, n, distribution);
        R.pop_back();
        completed.insert(i);
    }

    cout << "Number of cliques: " << count_cliques << endl;
    cout << "Size of largest maximal clique: " << max_size << endl;

    auto end = chrono::high_resolution_clock::now();
    cout << endl;
    cout << "Distribution of clique sizes:" << endl;
    cout << "Size Number of cliques" << endl;
    for (auto x: distribution) {
        cout << x.first << "    " << x.second << endl;
    }
    cout << endl;
    chrono::duration<double> duration = end - start;
    cout << "Time taken by algorithm: " << duration.count() << " seconds" << endl;


    return 0;
}
