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
uint n, m, key_cap, max_key, min_key;
uint *keys, *heads, *pres, *nexts;

uint *id, *degree, *core, *seq, max_core;
uint *pstart, *edges;


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


//time complexity: O(m)
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

            for(uint j = pstart[id]; j < pstart[id+1]; j++) {
                uint v = edges[j];
                if(keys[v] != n)
                    decrement(v, 1);
            }
        }
    }

    printf("Algorithm 1: peel:\n");
    printf("core(v) for each v （ V:\n");
    for(uint i = 0; i < n; i++)
        printf("core(%d) = %d\n", i, core[i]);
}


//time complexity: O(m)
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

        for(uint i = pstart[u]; i < pstart[u+1]; i++) {
            uint v = edges[i];
            if(d[v] != n) {
                d[v] = d[v]-1;
                if(d[v] == max_core-1)
                    Q.push(v);
            }
        }
    }
    printf("\nAlgorithm 2: k-core:\n");
    printf("the subgraph of G induced by vertices with d(.) >= %d\n", max_core);
    for(uint i = 0; i < n; i++) {
        if(d[i] != n) {
            printf("neighbor vertex of %d: ", i);
            for(uint j = pstart[i]; j < pstart[i+1]; j++) {
                uint v = edges[j];
                if(d[v] != n)
                    printf("%d ", v);
            }
            printf("\n");
        }
    }
}


uint *pre;

void UF_init() {
    pre = (uint *)malloc(sizeof(uint)*n);
    for(uint i = 0; i < n; i++)
        pre[i] = i;
}

uint UF_find(uint u) {
    return pre[u] == u ? u : pre[u] = UF_find(pre[u]);
}

uint UF_union(uint u, uint v) {
    uint fu = UF_find(u);
    uint fv = UF_find(v);

    if(fu != fv)
        pre[fu] = pre[fv] = min(fu, fv);
}

set<uint> ss[maxn];
set<uint>::iterator iter;
struct node {
    uint data;
    uint firstchild;
    uint nextsibling;
} tree[maxn];

void CoreHierarchy() {   //compute the core hierarchy tree of a graph
    UF_init();
    uint *vis = (uint *)malloc(sizeof(uint)*n);
    memset(vis, 0, sizeof(uint)*n);
    uint root;

    for(uint i = 0; i < n; i++) {
        ss[i].insert(i);
        tree[i].data = core[i];
        tree[i].firstchild = -1;
        tree[i].nextsibling = -1;
    }

    for(uint i = n; i > 0; i--) {
        uint u = seq[i-1];
        for(uint j = pstart[u]; j < pstart[u+1]; j++) {
            uint v = edges[j];
            if(vis[v]) {
                uint ru = UF_find(u);
                uint rv = UF_find(v);

                if(ru != rv) {
                    if(core[ru] == core[rv]) {
                        if(ru > rv)
                            swap(ru, rv);
                        ss[ru].insert(ss[rv].begin(), ss[rv].end());
                        pre[rv] = ru;
                        if(tree[rv].firstchild != -1) {
                            if(tree[ru].firstchild == -1)
                                tree[ru].firstchild = tree[rv].firstchild;
                            else {
                                uint p = tree[ru].firstchild;
                                while(tree[p].nextsibling != -1)
                                    p = tree[p].nextsibling;
                                tree[p].nextsibling = tree[rv].firstchild;
                            }
                        }

                        if(tree[rv].nextsibling != -1) {
                            uint p = ru;
                            while(tree[p].nextsibling != -1)
                                p = tree[p].nextsibling;
                            tree[p].nextsibling = tree[rv].nextsibling;
                        }
                    }
                    else {
                        if(tree[ru].firstchild == -1)
                            tree[ru].firstchild = rv;
                        else {
                            uint p = tree[ru].firstchild;
                            while(tree[p].nextsibling != -1)
                                p = tree[p].nextsibling;
                            tree[p].nextsibling = rv;
                        }
                    }
                    root = ru;
                }
            }
        }
        vis[u] = 1;
    }
    printf("\nAlgorithm 3: CoreHierarchy:\n");
    printf("A core hierarchy tree CoreHT of G:\n");
    free(vis);
}

vector<pair<pair<uint, uint>, uint> > CoreSPT;

void CoreSpanning() {    //compute a core spanning tree of a graph
    UF_init();
    uint *vis = (uint *)malloc(sizeof(uint)*n);
    memset(vis, 0, sizeof(uint)*n);

    for(uint i = n; i > 0; i--) {
        uint u = seq[i-1];
        for(uint j = pstart[u]; j < pstart[u+1]; j++) {
            uint v = edges[j];
            if(vis[v]) {
                uint ru = UF_find(u);
                uint rv = UF_find(v);

                if(ru != rv) {
                    UF_union(ru, rv);
                    CoreSPT.push_back(make_pair(make_pair(u, v), core[u]));
                }
            }
        }
        vis[u] = 1;
    }
    printf("\nAlgorithm 4: CoreSpanning:\n");
    printf("A core spanning tree CoreSPT of G\n");
    for(uint i = 0; i < CoreSPT.size(); i++) {
        printf("%d %d %d\n", CoreSPT[i].first.first, CoreSPT[i].first.second, CoreSPT[i].second);
    }
    free(vis);
}

uint Hindex(uint u, uint *core_up) {
    uint s = degree[u];
    uint *cnt = (uint *)malloc(sizeof(uint)*(s+1));
    memset(cnt, 0, sizeof(uint)*(s+1));
    for(uint i = pstart[u]; i < pstart[u+1]; i++) {
        uint v = edges[i];
        if(core_up[v] > s)
            cnt[s]++;
        else
            cnt[core_up[v]]++;
    }
    for(uint i = s; i > 0; i--) {
        if(cnt[i] >= i)
            return i;
        cnt[i-1] = cnt[i-1]+cnt[i];
    }
}


void CoreD_Local() {    //compute core numbers of vertices
    uint *core_up = (uint *)malloc(sizeof(uint)*n);
    for(uint i = 0; i < n; i++)
        core_up[i] = degree[i];
    uint update = 1;
    while(update) {
        update = 0;
        for(uint i = 0; i < n; i++) {
            uint old = core_up[i];
            core_up[i] = Hindex(i, core_up);
            if(core_up[i] != old)
                update = 1;
        }
    }
    printf("\nAlgorithm 5: CoreD-Local:\n");
    printf("core(v) of each vertex v （ V\n");
    for(uint i = 0; i < n; i++) {
        core[i] = core_up[i];
        printf("core(%d) = %d\n", i, core[i]);
    }
    free(core_up);
}

void CoreD_Local_opt() {    //compute core numbers of vertices
    uint *core_up = (uint *)malloc(sizeof(uint)*n);
    uint *cnt = (uint *)malloc(sizeof(uint)*n);
    uint *vis = (uint *)malloc(sizeof(uint)*n);
    memset(cnt, 0, sizeof(uint)*n);
    memset(vis, 0, sizeof(uint)*n);

    queue<uint> Q, Q1;
    for(uint i = 0; i < n; i++) {
        core_up[i] = degree[i];
        for(uint j = pstart[i]; j < pstart[i+1]; j++) {
            if(degree[edges[j]] >= degree[i])
                cnt[i]++;
        }
        if(cnt[i] < core_up[i]) {
            Q.push(i);
            vis[i] = 1;
        }
    }

    while(!Q.empty()) {
        while(!Q.empty()) {
            uint u = Q.front();
            Q.pop();
            vis[u] = 0;

            uint old = core_up[u];
            core_up[u] = Hindex(u, core_up);
            cnt[u] = 0;
            for(uint i = pstart[u]; i < pstart[u+1]; i++) {
                uint v = edges[i];
                if(core_up[v] >= core_up[u])
                    cnt[u]++;
                if(!vis[v] && core_up[u] < core_up[v] && core_up[v] <= old) {
                    if(cnt[v] == core_up[v])
                        Q1.push(v);
                    cnt[v]--;
                }
            }
        }
        while(!Q1.empty()) {
            uint x = Q1.front();
            Q1.pop();
            Q.push(x);
            vis[x] = 1;
        }
    }
    printf("\nAlgorithm 6: CoreD-Local-opt:\n");
    printf("core(v) of each vertex v （ V\n");
    for(uint i = 0; i < n; i++) {
        core[i] = core_up[i];
        printf("core(%d) = %d\n", i, core[i]);
    }
    free(core_up);
    free(cnt);
    free(vis);
}



void CoreD_IO() {     //I/O efficiently computed core numbers of vertices
    uint v_min = n-1, v_max = 0;
    uint *core_up = (uint *)malloc(sizeof(uint)*n);
    uint *cnt = (uint *)malloc(sizeof(uint)*n);
    memset(cnt, 0, sizeof(uint)*n);

    for(uint i = 0; i < n; i++) {
        core_up[i] = degree[i];
        for(uint j = pstart[i]; j < pstart[i+1]; j++) {
            if(degree[edges[j]] >= degree[i])
                cnt[i]++;
        }
        if(cnt[i] < core_up[i]) {
            if(i < v_min)
                v_min = i;
            v_max = i;
        }
    }

    while(v_min <= v_max) {
        uint v1_min = n-1, v1_max = 0;
        for(uint u = v_min; u <= v_max; u++) {
            if(core_up[u] <= cnt[u])
                continue;
            uint old = core_up[u];
            core_up[u] = Hindex(u, core_up);
            cnt[u] = 0;

            for(uint j = pstart[u]; j < pstart[u+1]; j++) {
                uint v = edges[j];
                if(core_up[v] >= core_up[u])
                    cnt[u]++;
                if(core_up[u] < core_up[v] && core_up[v] <= old) {
                    if(cnt[v] == core_up[v] && !(u < v && v <= v_max)) {
                        if(v < v1_min)
                            v1_min = v;
                        if(v > v1_max)
                            v1_max = v;
                    }
                    cnt[v]--;
                }
            }
        }
        v_min = v1_min;
        v_max = v1_max;
    }
    printf("\nAlgorithm 7: CoreD-IO:\n");
    printf("core(v) of each vertex v （ V\n");
    for(uint i = 0; i < n; i++) {
        core[i] = core_up[i];
        printf("core(%d) = %d\n", i, core[i]);
    }
    free(core_up);
    free(cnt);
}

void solve() {
    scanf("%d", &n);
    key_cap = n-1;

    heads = (uint *)malloc(sizeof(uint)*(key_cap+1));
    nexts = (uint *)malloc(sizeof(uint)*n);
    pres = (uint *)malloc(sizeof(uint)*n);
    keys = (uint *)malloc(sizeof(uint)*n);
    id = (uint *)malloc(sizeof(uint)*n);
    degree = (uint *)malloc(sizeof(uint)*n);
    pstart = (uint *)malloc(sizeof(uint)*(n+1));

    m = 0;
    for(uint i = 0; i < n; i++) {
        id[i] = i;
        scanf("%d", &degree[i]);
        pstart[i] = m;
        m += degree[i];
    }
    pstart[n] = m;

    edges = (uint *)malloc(sizeof(uint)*m);
    uint cnt = 0;
    for(uint i = 0; i < n; i++) {
        for(uint j = 0; j < degree[i]; j++)
            scanf("%d", &edges[cnt++]);
    }

    peel();
    k_core();
    CoreHierarchy();
    CoreSpanning();
    CoreD_Local();
    CoreD_Local_opt();
    CoreD_IO();

    free(heads);
    free(nexts);
    free(pres);
    free(keys);
    free(id);
    free(degree);
    free(core);
    free(seq);
    free(pstart);
    free(edges);
}

int main() {
    freopen("data.txt","r",stdin);
    freopen("result.txt","w",stdout);
    solve();
}
