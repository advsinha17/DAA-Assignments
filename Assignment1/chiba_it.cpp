#include <stdio.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <set>
#include <memory>
using namespace std;

long long count_cliques = 0;
long long max_size = 0;


struct State {
    int i;
    int val;
    unique_ptr<vector<int>> CMinusN;
    unique_ptr<vector<int>> NMinusC;
    unique_ptr<vector<int>> CintersectN;
};

void update(vector<int> &CHash, vector<int>& C, int v, int n, vector<set<int>>& adj, vector<int>& T, vector<int>& S, unordered_map<int, int>& mp, unordered_map<int, int>& distribution) {
    stack<State> st;

    State initialState;
    initialState.i = v;
    initialState.val = 0;
    initialState.CMinusN = make_unique<vector<int>>();
    initialState.NMinusC = make_unique<vector<int>>();
    initialState.CintersectN = make_unique<vector<int>>();
    st.push(std::move(initialState));

    while (!st.empty()) {
        State s = std::move(st.top());
        int i = s.i;
        st.pop();

        if (i == n) {
            count_cliques++;
            max_size = max(max_size, (long long)C.size());
            distribution[(int)C.size()]++;
            // cout << "found clique " << count_cliques << endl;
            continue;
        }

        if (s.val == 0) {
            auto CMinusN = make_unique<vector<int>>();
            auto NMinusC = make_unique<vector<int>>();
            auto CintersectN = make_unique<vector<int>>();

            for (auto x: adj[i]) {
                if (CHash[x] == 0) NMinusC->push_back(x);
                else CintersectN->push_back(x);
            }

            for (auto x: C) {
                if (adj[i].find(x) == adj[i].end()) CMinusN->push_back(x);
            }

            int cmnSize = CMinusN->size();
            s.CMinusN = std::move(CMinusN);
            s.NMinusC = std::move(NMinusC);
            s.CintersectN = std::move(CintersectN);

            if (cmnSize != 0) {
                s.val = 1;
                st.push(std::move(s));
                State nextState;
                nextState.i = i + 1;
                nextState.val = 0;
                nextState.CMinusN = make_unique<vector<int>>();
                nextState.NMinusC = make_unique<vector<int>>();
                nextState.CintersectN = make_unique<vector<int>>();
                st.push(std::move(nextState));
                continue;
            }

            s.val = 1;
        }

        if (s.val == 1) {
            vector<int>& CMinusN = *s.CMinusN;
            vector<int>& NMinusC = *s.NMinusC;
            vector<int>& CintersectN = *s.CintersectN;

            int cmnSize = CMinusN.size();
            int nmcSize = NMinusC.size();
            int cinSize = CintersectN.size();

            for (auto x: CintersectN) {
                for (auto y: adj[x]) {
                    if (CHash[y] == 0 && y != i) T[y]++;
                }
            }

            for (auto x: CMinusN) {
                for (auto y: adj[x]) {
                    if (CHash[y] == 0) S[y]++;
                }
            }

            bool flag = true;
            for (auto y: NMinusC) {
                if (y < i && T[y] == cinSize) {
                    flag = false;
                    break;
                }
            }

            sort(CMinusN.begin(), CMinusN.end());
            int p = cmnSize;
            for (int k = 0; k < p; k++) {
                for (auto y: adj[CMinusN[k]]) {
                    if (CHash[y] == 0 && y < i && T[y] == cinSize) {
                        if (y >= CMinusN[k]) S[y]--;
                        else {
                            if ((y < CMinusN[k]) && (k == 0 || y >= CMinusN[k - 1]) && (S[y] + k) == p) {
                                flag = false;
                            }
                        }
                    }
                }
            }

            if (flag) {
                if (!CintersectN.empty()) {
                    for (int y = 0; y < i; y++) {
                        if (CHash[y] == 0 && T[y] == cinSize && S[y] == 0) {
                            if (cmnSize == 0 || CMinusN.back() < y) {
                                flag = false;
                                break;
                            }
                        }
                    }
                } else if (cmnSize == 0 || CMinusN.back() < i - 1) flag = false;
            }

            for (auto x: CintersectN) {
                for (auto y: adj[x]) {
                    if (CHash[y] == 0 && y != i) {
                        T[y] = 0;
                    }
                }
            }

            for (auto x: CMinusN) {
                for (auto y: adj[x]) {
                    if (CHash[y] == 0) {
                        S[y] = 0;
                    }
                }
            }

            if (flag) {
                s.val = 2;
                for (auto x: CMinusN) {
                    if (CHash[x] == 1) {
                        CHash[x] = 0;
                        auto it = find(C.begin(), C.end(), x);
                        if (it != C.end()) C.erase(it);
                    }
                }

                if (CHash[i] == 0) {
                    CHash[i] = 1;
                    C.push_back(i);
                }

                st.push(std::move(s));
                State nextState;
                nextState.i = i + 1;
                nextState.val = 0;
                nextState.CMinusN = make_unique<vector<int>>();
                nextState.NMinusC = make_unique<vector<int>>();
                nextState.CintersectN = make_unique<vector<int>>();
                st.push(std::move(nextState));
                continue;
            }
        }

        if (s.val == 2) {
            if (CHash[i]) {
                CHash[i] = 0;
                auto it = find(C.begin(), C.end(), i);
                if (it != C.end()) C.erase(it);
            }

            vector<int>& CMinusN = *s.CMinusN;
            for (auto x: CMinusN) {
                if (CHash[x] == 0) {
                    CHash[x] = 1;
                    C.push_back(x);
                }
            }
        }
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
    vector<vector<int>> adj;
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
        
        if (adj.size() < n) {
            adj.resize(n);
            degrees.resize(n, 0);
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
        int mappedU = vertexMap[u];
        int mappedV = vertexMap[v];

        adj[mappedU].push_back(mappedV);
        adj[mappedV].push_back(mappedU);

    }
    infile.close();

    for (int i = 0; i < n; i++) {
        degrees[i] = adj[i].size();
    }
    vector<pair<int, int> > toSort;
    for (int i = 0; i < n; i++) {
        toSort.push_back({degrees[i], i});
    }
    sort(toSort.begin(), toSort.end());
    unordered_map<int, int> mp;
    unordered_map<int, int> rev_mp;
    for (int i = 0; i < n; i++) {
        mp[i] = toSort[i].second;
        rev_mp[toSort[i].second] = i;
    }
    vector<set<int> > mappedAdj(n);
    for (int i = 0; i < n; i++) {
        int originalVertex = mp[i];
        for (int neighbor : adj[originalVertex]) {
            mappedAdj[i].insert(rev_mp[neighbor]);
        }
    }
    vector<int> S(n, 0);
    vector<int> T(n, 0);
    unordered_map<int, int> distribution;
    vector<int> C;
    int v = 0;
    while (toSort[v].first == 0) {
        v++;
        count_cliques++;
    }
    C.push_back(v);
    vector<int> CHash(n, 0);
    CHash[v] = 1;
    update(CHash, C, v + 1, n, mappedAdj, T, S, mp, distribution);
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