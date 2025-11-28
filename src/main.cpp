#include "btbf.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <random>
#include <unordered_set>

static void printKey(const btbf::BTBF::SymKey& k, const char* label){
    std::cout << label << ": 0x";
    for (auto b : k) {
        std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b);
    }
    std::cout << std::dec << "\n";
}

int main(){
    std::mt19937_64 rng{ std::random_device{}() };
    int ns[] = {5, 10};
    int ts[] = {5};
    int ells[] = {5, 10};

    for(int n : ns){
        for(int t : ts){
            for(int ell : ells){
                std::cout << "n=" << n << ", t=" << t << ", ell=" << ell << "\n";
                btbf::BTBF::init();

                auto keygenOut = btbf::BTBF::keygen(n, t, ell, 256);

                std::uniform_int_distribution<int> distEll(1, ell);
                int j = distEll(rng); // position to encapsulate to

                // Encapsulate
                auto [k_enc, ct] = btbf::BTBF::encaps(keygenOut.pk, j);

                std::vector<int> J;
                J.reserve(t);
                std::unordered_set<int> used;
                std::uniform_int_distribution<int> distN(1, n);
                while (static_cast<int>(J.size()) < t) {
                    int idx = distN(rng);
                    if (used.insert(idx).second) {
                        J.push_back(idx);
                    }
                }

                std::vector<btbf::BTBF::GT> shares;
                shares.reserve(J.size());
                for (int idx : J) {
                    const auto& party = keygenOut.parties[idx - 1];
                    if(idx%2 == 0){
                    const auto& sk_left_j = party.left[j - 1]; // left key (b=0) for position j
                    auto di = btbf::BTBF::decShare(sk_left_j, ct.c0);
                    shares.push_back(di);
                    } else {
                        const auto& sk_right_j = party.right[j - 1]; // right key (b=1) for position j
                        auto di = btbf::BTBF::decShare(sk_right_j, ct.c1);
                        shares.push_back(di);
                    }
                }

                auto k_dec = btbf::BTBF::combine(J, shares);

                printKey(k_enc, "Encapsulated key");
                printKey(k_dec, "Recovered   key");

                bool equal = std::equal(k_enc.begin(), k_enc.end(), k_dec.begin());
                std::cout << "Keys match? " << (equal ? "YES" : "NO") << "\n";
            }
        }
    }
    return 0;
}
