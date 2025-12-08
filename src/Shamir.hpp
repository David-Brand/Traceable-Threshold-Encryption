#pragma once

#include <vector>

template<typename T>
std::vector<T> share(T a, int n, int t){
    std::vector<T> coef(t, std::rand());
    coef[0] = a;
    std::vector<T> s(n);

    for (int i = 0; i < n; i++) {
        int idx = i + 1;    // shares at positions 1..n
        T x = idx;

        T val = coef[0];
        T base = x;

        for (int k = 1; k < t; k++) {
            T term;
            term = coef[k] * base;

            val = val + term;

            base = base * x;
        }
        s[i] = val; // s_{i+1}
    }
    return s;
}