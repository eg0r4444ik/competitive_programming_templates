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
// суффиксный массив + массив lcp

void count_sort(vector<int> &p, vector<int> &c){
    int n = p.size();

    vector<int> cnt(n);
    for(auto x : c){
        cnt[x]++;
    }

    vector<int> p_new(n);
    vector<int> pos(n);
    pos[0] = 0;
    for(int i = 1; i < n; i++){
        pos[i] = pos[i-1] + cnt[i-1];
    }
    for(auto x : p){
        int i = c[x];
        p_new[pos[i]] = x;
        pos[i]++;
    }
    p = p_new;
}

pair<vector<ll>, vector<ll>> suffarray(string s){
    int n = s.length();
    vector<int> p(n);
    vector<int> c(n);
    vector<int> lcp(n);

    {
        vector<pair<char, int>> a(n);
        for(int i = 0; i < n; i++) a[i] = {s[i], i};
        sort(a.begin(), a.end());

        for(int i = 0; i < n; i++) p[i] = a[i].second;
        c[p[0]] = 0;
        for(int i = 1; i < n; i++){
            if(a[i].first == a[i-1].first){
                c[p[i]] = c[p[i-1]];
            }else{
                c[p[i]] = c[p[i-1]] + 1;
            }
        }
    }

    int k = 0;
    while ((1 << k) < n) {
        
        for(int i = 0; i < n; i++){
            p[i] = (p[i] - (1 << k) + n)%n;
        }

        count_sort(p, c);

        vector<int> c_new(n);
        c_new[p[0]] = 0;
        for(int i = 1; i < n; i++){
            pair<int, int> prev = {c[p[i-1]], c[(p[i-1] + (1 << k)) % n]};
            pair<int, int> curr = {c[p[i]], c[(p[i] + (1 << k)) % n]};
            if(prev == curr){
                c_new[p[i]] = c_new[p[i-1]];
            }else{
                c_new[p[i]] = c_new[p[i-1]] + 1;
            }
        }
        c = c_new;
        k++;
    }

    k = 0;
    for(int i = 0; i < n - 1; i++){
        int pi = c[i];
        int j = pi - 1 < 0 ? 0 : p[pi - 1];
        while(s[i + k] == s[j + k]) k++;
        lcp[pi] = k;
        k = max(k - 1, 0ll);
    }

    return {p, lcp};
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// z-функция

vector<ll> zFunc(string s){
    int n = s.length();
    vector<ll> z(s.length());

    ll l = 0;
    ll r = 0;
    z[0] = 0;

    for(int i = 0; i < s.length(); i++){
        if(i < r){
            z[i] = min(r-i+1, z[i-l]);
        }
        while(i + z[i] < s.length() && s[i + z[i]] == s[z[i]]){
            z[i]++;
        }
        if(r < z[i] + i - 1){
            l = i;
            r = i + z[i] - 1;
        }
    }

    return z;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// префикс-функция

vector<ll> prefFunc (string s) {
	int n = s.length();
	vector<ll> pref(n);

	for (int i = 1; i < n; i++) {
		int j = pref[i-1];
		while (j > 0 && s[i] != s[j])
			j = pref[j-1];
		if (s[i] == s[j])  j++;
		pref[i] = j;
	}

	return pref;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// хэши

ll pow(ll n, ll p, ll m){
    if(p == 0){
        return 1;
    }
    if(p == 1){
        return n;
    }
    if(p%2 == 0){
        ll res = pow(n, p/2, m)%m;
        return (res*res)%m;
    }else{
        ll res = pow(n, p/2, m);
        return (((res*res)%m)*n)%m;
    }
}

vector<ll> hashs(string s){
    ll n = s.length();
    ll p = 31;
    ll mod = 1000000007;
 
    vector<ll> hashs(n);
    hashs[0] = s[0] - 'a' + 1;
    for(int i = 1; i < n; i++){
        hashs[i] = (hashs[i-1] + (s[i] - 'a' + 1)*pow(p, i, mod)%mod)%mod;
    }

    return hashs;
}

ll substrHash(int l, int r, vector<ll> hashs){
    if(l == 0){
        return hashs[r];
    }

    ll p = 31;
    ll mod = 1000000007;

    return ((hashs[r] - hashs[l-1] + mod)%mod*pow(pow(p, mod-2, mod), l, mod))%mod;
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// бор

struct Trie{

    struct Node{
        vector<int> to;
        bool term;

        Node() : to(26), term(false) {}
    };

    vector<Node> data;
    int root;
    int free;

    Trie() : data(10000000), root(0), free(1) {
        data[0] = Node();
    }

    void add(string s){
        int curr = root;
        for(char c : s) {
            if(data[curr].to[c-'a'] == 0) {
                data[curr].to[c-'a'] = free++;
                curr = data[curr].to[c-'a'];
                data[curr] = Node();
            } else {
                curr = data[curr].to[c-'a'];
            }
        }
        data[curr].term = true;
    }

    bool find(string s){
        int curr = root;
        for(char c : s) {
            if(data[curr].to[c-'a'] == 0) {
                return false;
            }
            curr = data[curr].to[c-'a'];
        }
        if(data[curr].term){
            return true;
        } else{
            return false;
        }
    }
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Ахо-Корасик