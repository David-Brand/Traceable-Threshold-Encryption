// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: GPL-3.0-or-later

#include "btbf.hpp"
#include "TTT.hpp"
#include "test_functions.cpp"

#include <chrono>

int main(){
    using namespace ttt;

    bool testBTDDH = false;
    bool testBTBF = true;

    std::mt19937_64 rng{ std::random_device{}() };
    int ns[] = {100, 1000};
    int ts[] = {25, 100};

    int lambdas[] = {128};
    double es[] = {0.99, 0.5};

    for(int n : ns){
        for(int t : ts){
            for(int lambda : lambdas){
                for(double e : es){
                    if(testBTDDH){
                        testTTT_BTDDH(n, t, lambda, e);
                    }
                    if(testBTBF){
                        testTTT_BTBF(n, t, lambda, e);
                    }
                }
            }
        }
    }
    return 0;
}
