#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <set>
#include <chrono>
using namespace std;


long long count_cliques = 0;
long long max_size = 0;

void expand(unordered_set<int> &subg, unordered_set<int> &cand,vector<set<int> > &adj, vector<int>& Q, int n, unordered_map<int, int> &distribution){
    if(subg.empty()){
        distribution[(int)Q.size()]++;
        count_cliques++;
        max_size = max(max_size, (long long)Q.size());
        return;
    }
    int pivot = -1;
    int pivotSize = -1;
    for (auto i: subg) {
        int count = 0;
        if (adj[i].size() < cand.size()) {
            for (auto y: adj[i]) {
                if (cand.find(y) != cand.end()) count++;
            }
        } else {
            for (auto y: cand) {
                if (adj[i].find(y) != adj[i].end()) count++;
            }
        }

        if (count > pivotSize) {
            pivot = i;
            pivotSize = count;
        }
    }
    for (auto it = cand.begin(); it != cand.end(); ) {
        int x = *it;
        if (adj[pivot].find(x) != adj[pivot].end()) {
            ++it;
            continue;
        } 
        unordered_set<int> newSubg;
        unordered_set<int> newCand;
        if (cand.size() < adj[x].size()) {
            for (auto y: cand) {
                if (adj[x].find(y) != adj[x].end()) newCand.insert(y);
            }
        } else {
            for (auto y: adj[x]) {
                if (cand.find(y) != cand.end()) newCand.insert(y);
            }
        }
        
        if (subg.size() < adj[x].size()) {
            for (auto y: subg) {
                if (adj[x].find(y) != adj[x].end()) newSubg.insert(y);
            }
        } else {
            for (auto y: adj[x]) {
                if (subg.find(y) != subg.end()) newSubg.insert(y);
            }
        }
        Q.push_back(x);
        expand(newSubg, newCand, adj, Q, n, distribution);
        Q.pop_back();
        it = cand.erase(it);
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
        }

        int mappedU = vertexMap[u];
        int mappedV = vertexMap[v];

        adj[mappedU].insert(mappedV);
        adj[mappedV].insert(mappedU);
    }

    infile.close();

    unordered_set<int> subg(n);
    unordered_set<int> cand(n);
    unordered_map<int, int> distribution;
    for (int i = 0; i < n; i++) {;
        subg.insert(i);
        cand.insert(i);
    }
    vector<int> Q;
    expand(subg, cand, adj, Q, n, distribution);

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
