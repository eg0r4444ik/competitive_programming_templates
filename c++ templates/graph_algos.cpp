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
// lca

struct LCA{

    struct Node{
        int val;
        bool isVisited;
        vector<Node> neighbours;

        Node(ll val): val(val), isVisited(false){}
    };

    vector<vector<int>> up;
    vector<int> dist;
    vector<Node> tree;

    LCA(int n) {
        up.resize(n+1, vector<int>(31));
        dist.resize(n+1);
        tree.resize(n+1);

        dist[1] = 0;
        for(int i = 0; i <= 30; i++){
            up[1][i] = 0;
        }

        for(int i = 0; i < n+1; i++){
            tree[i] = Node(i);
        }
    }

    void add(int a, int b){
        tree[a].neighbours.push_back(tree[b]);
        tree[b].neighbours.push_back(tree[a]);

        dist[b] = dist[a]+1;
        up[b][0] = a;
        for(int j = 1; j <= 30; j++){
            up[b][j] = up[up[b][j-1]][j-1];
        }
    }

    int lca(int a, int b){
        if (dist[a] < dist[b]){
            int tmp = a;
            a = b;
            b = tmp;
        }

        for (int j = 30; j >= 0; j--){
            if (dist[up[a][j]] >= dist[b] && up[a][j] != 0){
                a = up[a][j];
            }
        }
        if(a == b){
            return b;
        }else {
            for (int j = 30; j >= 0; j--) {
                if (up[a][j] == up[b][j]) {
                    continue;
                } else {
                    a = up[a][j];
                    b = up[b][j];
                }
            }
            return up[a][0];
        }
    }
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// паросочетания

vector<int> mt;
int n, k;
vector<bool> used;
vector<Node> graph;

struct Node{
    int val;
    vector<Node> neighbours;

    Node(int val): val(val) {}
};

void solve() {
    // инициализация графа
    cin >> n >> k;
    graph.resize(n+k);
    used.resize(n+k);
    mt.resize(k);
    for(int i = 0; i < n+k; i++){
        graph[i] = Node(i);
    }
    for(int i = 0; i < n; i++){
        int num;
        cin >> num;
        while(num != 0){
            graph[i].neighbours.add(graph[n+num-1]);
            graph[n+num-1].neighbours.add(graph[i]);
            cin >> num;
        }
    }

    // поиск максимального паросочетания
    fill(mt.begin(), mt.end(), -1);
    fill(used.begin(), used.end(), false);


    for (int i = 0; i < n; i++) {
        if(tryKuhn(graph[i])){
            fill(used.begin(), used.end(), false);
        }
    }

    vector<string> res;
    for(int i = 0; i < mt.size(); i++){
        if(mt[i] != -1){
            res.push_back((mt[i]+1) + " " + (i+1));
        }
    }

    cout << res.size(); // размер максимального паросочетания
    for(string s : res){
        cout << s << nl; // пара вершин в паросочетании
    }
}

bool tryKuhn(Node curr){
    if(used[curr.val]){
        return false;
    }
    used[curr.val] = true;
    for(Node node : curr.neighbours){
        int idx = node.val - n;
        if(mt[idx] == -1){
            mt[idx] = curr.val;
            return true;
        }else if(tryKuhn(graph[mt[idx]])){
            mt[idx] = curr.val;
            return true;
        }
    }

    return false;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// потоки

struct FlowGraph {
    struct Edge {
        ll to;
        ll flow, cost, cap;

        Edge(ll to, ll cap, ll flow, ll cost) : to(to), cap(cap), flow(flow), cost(cost) {}
    };

    struct Pair {
        long flow, cost;

        Pair(ll flow, ll cost) : flow(flow), cost(cost) {}
    };

    vector<Edge> edges;
    vector<vector<ll>> graph;
    vector<bool> used;
    vector<ll> lst;

    vector<ll> p;
    vector<ll> pe;
    vector<ll> dist;
    int s, t, n;

    void init(int k) {
        n = k;
        graph.resize(n);
        used.resize(n);
        lst.resize(n);
        p.resize(n);
        pe.resize(n);
        dist.resize(n);
        s = 0;
        t = n - 1;
    }

    void addEdge(int from, int to, ll cap, ll flow, ll cost) {
        graph[from].push_back(edges.size());
        edges.push_back(Edge(to, cap, flow, cost));
        graph[to].push_back(edges.size());
        edges.push_back(Edge(from, 0, -flow, -cost));
    }

    void addBoundedEdge(int from, int to, ll l, ll r, ll flow, ll cost) {
        graph[from].push_back(edges.size());
        edges.push_back(Edge(to, r - l, flow, cost));
        graph[to].push_back(edges.size());
        edges.push_back(Edge(from, 0, -flow, -cost));

        graph[from].push_back(edges.size());
        edges.push_back(Edge(n - 1, l, flow, cost));
        graph[n - 1].push_back(edges.size());
        edges.push_back(Edge(from, 0, -flow, -cost));

        graph[0].push_back(edges.size());
        edges.push_back(Edge(to, l, flow, cost));
        graph[to].push_back(edges.size());
        edges.push_back(Edge(0, 0, -flow, -cost));
    }

    ll res(int x) {
        return edges[x].cap - edges[x].flow;
    }

    ll dfs(int x, ll f, int k) {
        if (used[x]) { return 0; }
        if (x == t) { return f; }
        used[x] = true;

        for (int e : graph[x]) {
            if (res(e) < k) { continue; }
            ll pushed = dfs(edges[e].to, min(f, res(e)), k);
            if (pushed != 0) {
                edges[e].flow += pushed;
                edges[e ^ 1].flow -= pushed;
                return pushed;
            }
        }
        return 0;
    }

    Pair augment() {
        spfa();
        if (p[t] == -1) return Pair(0, 0);
        ll mf = 1e18;
        int mc = 0;
        int curr = t;
        while (curr != s) {
            int e = pe[curr];
            mf = min(mf, res(pe[curr]));
            mc += edges[e].cost;
            curr = p[curr];
        }
        curr = t;
        while (curr != s) {
            int e = pe[curr];
            edges[e].flow += mf;
            edges[e ^ 1].flow -= mf;
            curr = p[curr];
        }
        return Pair(mf, mc * mf);
    }

    void spfa() {
        fill(p.begin(), p.end(), -1);
        fill(pe.begin(), pe.end(), -1);
        fill(dist.begin(), dist.end(), 1e9);
        for (int i = 0; i < n; i++) {
            p[i] = -1;
            pe[i] = -1;
            dist[i] = 1e9;
        }
        dist[s] = 0;
        queue<ll> qq;
        qq.push(s);
        vector<bool> inq (n);
        inq[s] = true;
        while (!qq.empty()) {
            int k = qq.front();
            qq.pop();
            inq[k] = false;
            for (int e : graph[k]) {
                if (res(e) == 0) continue;
                int to = edges[e].to;
                ll w = edges[e].cost;
                if (dist[to] > dist[k] + w) {
                    dist[to] = dist[k] + w;
                    p[to] = k;
                    pe[to] = e;
                    if (!inq[to]) {
                        qq.push(to);
                        inq[to] = true;
                    }
                }
            }
        }
    }

    bool bfs() {
        for (int i = 0; i < n; i++) {
            dist[i] = 1e9;
        }
        dist[s] = 0;
        queue<ll> qq;
        qq.push(s);
        while (!qq.empty()) {
            int k = qq.front();
            qq.pop();
            for (int e : graph[k]) {
                if (res(e) == 0) continue;
                int to = edges[e].to;
                if (dist[to] > dist[k] + 1) {
                    dist[to] = dist[k] + 1;
                    qq.push(to);
                }
            }
        }

        if (dist[t] != 1e9) {
            return true;
        }
        return false;
    }

    ll dinicDfs(int x, long mx) {
        if (x == t) return mx;
        int sum = 0;
        for (; lst[x] < graph[x].size(); lst[x]++) {
            int e = graph[x][lst[x]];
            int to = edges[e].to;
            ll r = res(e);
            if (r == 0 || dist[to] != dist[x] + 1) continue;
            ll pushed = dinicDfs(to, min(mx - sum, r));
            sum += pushed;
            edges[e].flow += pushed;
            edges[e ^ 1].flow -= pushed;
            if (sum == mx) break;
        }
        return sum;
    }

    ll dinic() {
        ll flow = 0;
        while (true) {
            if (!bfs()) break;
            fill(lst.begin(), lst.end(), 0);
            while (true) {
                long f = dinicDfs(s, 1e9);
                if (f == 0) break;
                flow += f;
            }
        }
        return flow;
    }

    ll maxFlow() {
        ll ans = 0;
        int k = 1 << 30;
        while (k >= 1) {
            while (true) {
                fill(used.begin(), used.end(), false);
                ll f = dfs(s, 1e9, k);
                if (f == 0) {
                    break;
                }
                ans += f;
            }
            k /= 2;
        }
        return ans;
    }

    Pair minCostMaxFlow() {
        ll flow = 0;
        ll cost = 0;
        while (true) {
            Pair pair = augment();
            if (pair.flow == 0) {
                break;
            }
            flow += pair.flow;
            cost += pair.cost;
        }
        return Pair(flow, cost);
    }
};