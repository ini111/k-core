#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>
#include <cmath>
#include <algorithm>
#include <queue>
#include <set>
#include <map>
#include <vector>
#include <iostream>
using namespace std;
typedef unsigned int uint;
const uint maxn = 1e3+10;

//the space complexity of ListLinearHeap is 3*n+key_cap+O(1)
uint n, key_cap, max_key, min_key;
uint *keys, *heads, *pres, *nexts;

void insert(uint id, uint key) {
    keys[id] = key, pres[id] = n, nexts[id] = heads[key];
    if(heads[key] != n)
        pres[heads[key]] = id;
    heads[key] = id;

    if(key < min_key) min_key = key;
    if(key > max_key) max_key = key;
}

//time complexity: O(n+key_cap)
void init(uint _n, uint _key_cap, uint *_ids, uint *_keys) {
    max_key = 0, min_key = _key_cap;
    for(int i = 0; i <= _key_cap; i++)
        heads[i] = n;

    for(int i = 0; i <= _n-1; i++)
        insert(_ids[i], _keys[i]);
}

uint remove(uint id) {
    if(pres[id] == n) {
        heads[keys[id]] = nexts[id];
        if(nexts[id] != n)
            pres[nexts[id]] = n;
    }
    else {
        nexts[pres[id]] = nexts[id];
        if(nexts[id] != n)
            pres[nexts[id]] = pres[id];
    }
    return keys[id];
}

//best time complexity: O(1), worst time complexity: O(max_key)
bool pop_min(uint &id, uint &key) {
    while(min_key <= max_key && heads[min_key] == n)
        min_key = min_key+1;
    if(min_key > max_key)
        return false;
    else {
        id = heads[min_key], key = min_key;
        remove(id);
        return true;
    }
}

uint decrement(uint id, uint dec) {
    uint key = remove(id);
    key = key-dec;
    insert(id, key);
    return key;
}

uint get_key(uint id) {
    return keys[id];
}

uint increment(uint id, uint inc) {

}

bool get_max(uint &id, uint &key) {

}

bool pop_max(uint &id, uint &key) {

}

bool get_min(uint &id, uint &key) {

}


vector<uint> adj[maxn];
uint *id, *degree, *core, *seq, max_core;

//time complexity: m+O(n)
void peel() {      //computer core number of all vertices
    init(n, key_cap, id, degree);
    core = (uint *)malloc(sizeof(uint)*n);
    seq = (uint *)malloc(sizeof(uint)*n);

    max_core = 0;
    for(uint i = 0; i < n; i++) {
        uint id, key;
        if(pop_min(id, key)) {
            seq[i] = id;
            if(keys[id] > max_core)
                max_core = keys[id];
            core[id] = max_core;
            keys[id] = n;

            for(uint j = 0; j < adj[id].size(); j++) {
                uint v = adj[id][j];
                if(keys[v] != n)
                    decrement(v, 1);
            }
        }
    }

    printf("core(v) for each vector v ¡Ê V:\n");
    for(uint i = 0; i < n; i++)
        printf("core(%d) = %d\n", i, core[i]);
}


//time complexity: O(n+m)
void k_core() {            //computer k-core
    uint *d = (uint *)malloc(sizeof(uint)*n);
    queue<uint> Q;

    for(uint i = 0; i < n; i++) {
        d[i] = degree[i];
        if(d[i] < max_core) {
            Q.push(i);
        }
    }

    while(!Q.empty()) {
        uint u = Q.front();
        Q.pop();
        d[u] = n;

        for(uint i = 0; i < adj[u].size(); i++) {
            uint v = adj[u][i];
            if(d[v] != n) {
                d[v] = d[v]-1;
                if(d[v] == max_core-1)
                    Q.push(v);
            }
        }
    }
    printf("\nthe subgraph of G induced by vertices with d(.) >= max_core\n");
    for(uint i = 0; i < n; i++) {
        if(d[i] != n) {
            printf("neighbor vertex of %d: ", i);
            for(uint j = 0; j < adj[i].size(); j++) {
                uint v = adj[i][j];
                if(d[v] != n)
                    printf("%d ", v);
            }
            printf("\n");
        }
    }
}


////the space complexity of the disjoint-set data structure is 2*n+O(1)
//uint *parent, *order;
//
//void UF_init(uint _n) {
//    for(uint i = 0; i < _n; i++) {
//        parent[i] = i;
//        order[i] = 0;
//    }
//}
//
//void UF_add(uint u, uint v) {
//    parent[u] = v;
//    order[u] = 0;
//}
//
//uint UF_find(uint u) {
//    uint res = u;
//    while(parent[id] != res)
//        res = parent[res];
//    while(parent[u] != res) {
//        uint tmp = parent[u];
//        parent[u] = res;
//        u = tmp;
//    }
//    return res;
//}
//
//uint UF_union(uint u, uint v) {
//    uint fu = UF_find(u);
//    uint fv = UF_find(v);
//
//    if(fu == fv)
//        return fu;
//
//    uint res;
//    if(order[fu] > order[fv]) {
//        res = fu;
//        parent[fv] = fu;
//    }
//    else {
//        res = fv;
//        parent[fu] = fv;
//        if(order[fu] == order[fv])
//            order[fv]++;
//    }
//    return res;
//}
//
//void CoreHierarchy() {   //compute the core hierarchy tree of a graph
//
//}

void solve() {
    scanf("%d", &n);
    key_cap = n-1;

    heads = (uint *)malloc(sizeof(uint)*(key_cap+1));
    nexts = (uint *)malloc(sizeof(uint)*n);
    pres = (uint *)malloc(sizeof(uint)*n);
    keys = (uint *)malloc(sizeof(uint)*n);
    id = (uint *)malloc(sizeof(uint)*n);
    degree = (uint *)malloc(sizeof(uint)*n);

    for(uint i = 0; i < n; i++) {
        id[i] = i;
        scanf("%d", &degree[i]);
    }

    uint temp;
    for(uint i = 0; i < n; i++) {
        for(uint j = 0; j < degree[i]; j++) {
            scanf("%d", &temp);
            adj[i].push_back(temp);
        }
    }

    peel();
    k_core();

    free(heads);
    free(nexts);
    free(pres);
    free(keys);
    free(id);
    free(degree);
    free(core);
    free(seq);
}

int main() {
    freopen("data.txt","r",stdin);
    freopen("result.txt","w",stdout);
    solve();
}
