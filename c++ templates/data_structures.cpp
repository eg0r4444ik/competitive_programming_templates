#include<iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <complex>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <queue>
#include <stack>
#include <string>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#define int long long
#define ll long long
#define nl "\n"

using namespace std;

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// дерево отрезков

struct SegmentTree{

    SegmentTree(){}

    struct Node{
        ll modify;
        ll function;

        Node(){}

        Node(ll modify, ll function): modify(modify), function(function){}
    };

    ll NEUTRAL_ELEMENT = 1e18;
    ll NO_OPERATION = -1e18;
    ll MOD = 1000000007;

    vector<Node> tree;
    int size;

    ll op_modify(ll a, ll b, int len){
        if(b == NO_OPERATION){
            return a;
        }
        return b;
    }

    ll op_func(ll a, ll b){
        return min(a, b);
    }

    void propagate(int x, int lx, int rx){
        if(tree[x].modify == NO_OPERATION || rx-lx == 1) return;
        int m = (lx+rx)/2;
        tree[2*x+1].modify = op_modify(tree[2*x+1].modify, tree[x].modify, 1);
        tree[2*x+1].function = op_modify(tree[2*x+1].function, tree[x].modify, m-lx);
        tree[2*x+2].modify = op_modify(tree[2*x+2].modify, tree[x].modify, 1);
        tree[2*x+2].function = op_modify(tree[2*x+2].function, tree[x].modify, rx-m);
        tree[x].modify = NO_OPERATION;
    }

    void init(int n){
        size = 1;
        while(size < n){
            size *= 2;
        }

        tree.resize(2*size-1);
    }

    void build(vector<ll> a, int x, int lx, int rx){
        if(rx-lx == 1){
            if(lx < a.size()) {
                tree[x] = Node(NO_OPERATION, a[lx]);
            }else{
                tree[x] = Node(NO_OPERATION, 0);
            }
        }else{
            int m = (lx+rx)/2;
            build(a, 2*x+1, lx, m);
            build(a, 2*x+2, m, rx);
            tree[x] = Node(NO_OPERATION, op_func(tree[2*x+1].function, tree[2*x+2].function));
        }
    }

    void build(vector<ll> a){
        init(a.size());
        build(a, 0, 0, size);
    }


    void modify(int l, int r, ll v, int x, int lx, int rx){
        propagate(x, lx, rx);
        if(l >= rx || r <= lx){
            return;
        }
        if(lx >= l && rx <= r){
            tree[x].modify = op_modify(tree[x].modify, v, 1);
            tree[x].function = op_modify(tree[x].function, v, rx-lx);
            return;
        }
        int m = (lx+rx)/2;
        modify(l, r, v, 2*x+1, lx, m);
        modify(l, r, v, 2*x+2, m, rx);
        tree[x].function = op_func(tree[2*x+1].function, tree[2*x+2].function);
        tree[x].function = op_modify(tree[x].function, tree[x].modify, rx-lx);
    }

    void modify(int l, int r, ll v){
        modify(l, r, v, 0, 0, size);
    }

    ll func(int l, int r, int x, int lx, int rx){
        propagate(x, lx, rx);
        if(l >= rx || r <= lx){
            return NEUTRAL_ELEMENT;
        }
        if(lx >= l && rx <= r){
            return tree[x].function;
        }
        int m = (lx+rx)/2;
        ll res1 = func(l, r, 2*x+1, lx, m);
        ll res2 = func(l, r, 2*x+2, m, rx);
        return op_func(res1, res2);
    }

    ll func(int l, int r){
        return func(l, r, 0, 0, size);
    }
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// разреженная таблица

struct SparseTable{
    vector<vector<ll>> table;
    vector<int> pows;

    SparseTable(vector<ll> a) {
        int n = a.size();
        table.resize(n, vector<ll>(20));
        for(int i = 0; i < n; i++){
            table[i][0] = a[i];
        }

        for(int k = 1; k < 20; k++){
            for(int i = 0; i + (1 << k) - 1 < n; i++){
                table[i][k] = min(table[i][k-1], table[i + (1 << (k-1))][k-1]);
            }
        }

        pows.resize(n+1);
        int curr = 0;
        for(int i = 0; i <= n; i++){
            while(1 << (curr+1) <= i){
                curr++;
            }
            pows[i] = curr;
        }
    }

    ll query(int l, int r){
        int k = pows[r-l+1];
        return min(table[l][k], table[r - (1 << k) + 1][k]);
    }
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// дерево Фенвика (на сумму)

struct FenvikTree {

    vector<ll> fenv;
    
    ll fenvSum(int idx){
        int start = idx&(idx+1);
        long res = fenv[idx];
        if(start != 0){
            res += fenvSum(start-1);
        }
        return res;
    }

    void init(vector<ll> a){
        int n = a.size();
        vector<ll> pref(n);
        pref[0] = a[0];
        for(int i = 1; i < n; i++){
            pref[i] = pref[i-1] + a[i];
        }

        fenv.resize(n);
        fenv[0] = pref[0];
        for(int i = 1; i < n; i++){
            fenv[i] = (i & (i + 1)) == 0 ? pref[i] : pref[i] - pref[(i & (i + 1)) - 1];
        }
    }
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// система непересекающихся множеств

struct DSU{
    vector<int> p;
    vector<int> r;

    DSU(int n) {
        p.resize(n+1);
        r.resize(n+1);
        for(int i = 0; i < n+1; i++){
            p[i] = i;
            r[i] = 1;
        }
    }

    int find(int x){
        if(x == p[x]){
            return x;
        }
        p[x] = find(p[x]);
        return p[x];
    }

    void un(int x, int y) {
        x = find(x);
        y = find(y);
        if (r[x] > r[y]) {
            int t = x;
            x = y;
            y = t;
        }
        if(x != y) {
            r[y] += r[x];
            p[x] = y;
        }
    }

    bool get(int x, int y){
        return find(x) == find(y);
    }
};