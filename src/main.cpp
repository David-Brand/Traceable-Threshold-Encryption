#include "btbf.hpp"
#include "TTT.hpp"

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

static void prntLength(int N, int c, int lambda, double del){
    float delta = (1/2 - del) / (1/2 - 2/std::sqrt(lambda));
    double e = 1 / std::pow(2, lambda);
    int n_ = 2 * c;
    double valL = 2.0 * c * std::log((2.0 * N) / e);
    std::size_t L_ = static_cast<std::size_t>(std::ceil(valL));
    double valD = ((4 * std::pow(n_, 2))/std::pow(1-delta, 2)) * std::log((4.0 * n_ * L_) / e);
    std::size_t d_ = static_cast<std::size_t>(std::ceil(valD));
    auto tole = L_ * d_ * (n_-1);
    std::cout << "fingerprinting length: " << tole << "\n";
}

int main(){
    using namespace ttt;

    std::mt19937_64 rng{ std::random_device{}() };
    int ns[] = {10};
    int ts[] = {5};

    int lambdas[] = {16};
    double es[] = {0.1};

    int ells[] = {5, 10};

    for(int n : ns){
        for(int t : ts){
            for(int lambda : lambdas){
                for(double e : es){
                    TTT E;
                    auto keygenOut = E.keygen(n, t, lambda, e);

                    std::cout << "coalition: [";
                    std::vector<int> J;
                    J.reserve(t);
                    std::unordered_set<int> used;
                    std::uniform_int_distribution<int> distN(0, n-1);
                    while (static_cast<int>(J.size()) < t) {
                        int idx = distN(rng);
                        if (used.insert(idx).second) {
                            J.push_back(idx);
                            std::cout << idx << ", ";
                        }
                    }
                    std::cout << "]\n";

                    auto D = [&J, &keygenOut](const btbf::BTBF::Ciphertext& c, const btbf::BTBF::SymKey& k) {
                        // Simulated decryption oracle
                        std::vector<btbf::BTBF::GT> shares;
                        shares.reserve(J.size());
                        for (int idx : J) {
                            shares.push_back(TTT::dec(keygenOut.parties[idx][c.j], c));
                        }
                        auto k_rec = TTT::combine(J, shares);
                        return std::equal(k.begin(), k.end(), k_rec.begin());
                    };

                    bool inC = false;
                    std::vector<std::size_t> accused = E.trace(keygenOut.pk, keygenOut.tk, D);
                    std::cout << "accused: [";
                    for (std::size_t i = 0; i < accused.size(); ++i) {
                        std::cout << accused[i] << (i + 1 < accused.size() ? ", " : "");
                        for(int j : J){
                            if(static_cast<std::size_t>(j) == accused[i]){
                                inC = true;
                            }
                        }
                    }
                    std::cout << "] " << (inC ? " (in coalition)" : " (not in coalition)") << "\n";
                }
            }

            /*
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
            }*/
        }
    }
    return 0;
}
