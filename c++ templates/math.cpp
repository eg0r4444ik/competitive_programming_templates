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
#define cd complex<double>
#define nl "\n"

using namespace std;

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// fft

vector<cd> fft(vector<cd> a, bool invert){
    int n = (int) a.size();
    if(n == 1){
        return a;
    }

    vector<cd> a0 (n/2);
    vector<cd> a1 (n/2);

    for(int i = 0; i < n/2; i++){
        a0[i] = a[2*i];
        a1[i] = a[2*i+1];
    }

    vector<cd> y0 = fft(a0, invert);
    vector<cd> y1 = fft(a1, invert);

    double ang = 2*M_PI/n * (invert ? -1 : 1);
    cd w (1);
    cd wn (cos(ang), sin(ang));
    for (int i = 0; i < n/2; i++) {
        a[i] = y0[i] + w*y1[i];
        a[i+n/2] = y0[i] - w*y1[i];
        if (invert)
            a[i] /= 2, a[i+n/2] /= 2;
        w *= wn;
    }

    return a;
}

// перемножение двух полиномов
vector<ll> multiply(vector<ll> a, vector<ll> b){
    vector<cd> fa (a.begin(), a.end()), fb (b.begin(), b.end());
    int n = 1;
    while (n < max (a.size(), b.size())) n <<= 1;
    n <<= 1;

    fa.resize (n), fb.resize (n);
    fa = fft (fa, false);
    fb = fft (fb, false);

    for (int i = 0; i < n; i++) {
        fa[i] = fa[i]*fb[i];
    }
    fa = fft (fa, true);

    vector<ll> res (n);
    for (int i = 0; i < n; i++) {
        res[i] =  (ll) round(fa[i].real());
    }

    return res;
}

// возведение полинома в степень
vector<ll> polynomialPow(vector<ll> a, int pow){
    if(pow == 1){
        return a;
    }
    if(pow%2 == 0){
        auto res = polynomialPow(a, pow/2);
        return multiply(res, res);
    }else{
        auto res = polynomialPow(a, pow/2);
        return multiply(a, multiply(res, res));
    }
}

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Расширенный алгоритм Евклида

struct Gcd {
    int x, y, g;

    Gcd(int x, int y, int g): x(x), y(y), g(g) {}

    Gcd gcd(int a, int b) {
        if (b == 0) {
            return Gcd(1, 0, a);
        }
        Gcd g = gcd(b, a % b);
        int x = g.y;
        int y = g.x - g.y * (a / b);
        return Gcd(x, y, g.g);
    }

    // X0 = Xg*(c/g)
    // Y0 = Yg*(c/g)
    // X = X0 + k*(b/g)
    // Y = Y0 - k*(a/g)
};

//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Вероятностный тест Ферма на простоту числа

ll mul(ll a, ll b, ll mod) {
    if (b == 1)
        return a;
    if (b % 2 == 0) {
        long long t = mul(a, b / 2, mod);
        return (2 * t) % mod;
    }
    return (mul(a, b - 1, mod) + a) % mod;
}
 
ll pows(ll num, ll pow, ll mod) {
    if (pow == 0)
        return 1;
    if (pow % 2 == 0) {
        long long t = pows(num, pow / 2, mod);
        return mul(t, t, mod) % mod;
    }
    return (mul(pows(num, pow - 1, mod), num, mod)) % mod;
}
 
ll gcd(ll a, ll b) {
    if (b == 0)
        return a;
    return gcd(b, a % b);
}
 
bool ferma(ll x) {
    if (x == 2)
        return true;
    srand(time(NULL));
    for (int i = 0; i < 10; i++) {
        ll a = (rand() % (x - 2)) + 2;
        if (gcd(a, x) != 1)
            return false;
        ll res = pows(a, x - 1, x);
        if (res != 1)
            return false;
    }
    return true;
}
