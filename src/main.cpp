#include "btbf.hpp"
#include "TTT.hpp"
#include "BTIB.hpp"
#include "CCA_BT_KEM.hpp"

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

void testBTIB(int n, int t, int ell){
    using namespace btib;
    std::mt19937_64 rng{ std::random_device{}() };
    std::cout << "n=" << n << ", t=" << t << ", ell=" << ell << "\n";
    BTIB_KEM1 b = BTIB_KEM1();
    Fr a;
    a.setByCSPRNG();
    G1 ID;
    G1::mul(ID, b.g1_, a);

    auto keygenOut = b.keygen(n, t, ell, 128);

    std::uniform_int_distribution<int> distEll(0, ell - 1);
    int j = distEll(rng);

    auto [k_enc, ct] = b.enc(keygenOut.pk, j, ID);

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

    std::vector<std::pair<bool, std::shared_ptr<btib::idkShare>>> shares;
    shares.reserve(J.size());
    for (int idx : J) {
        const auto& party = keygenOut.parties[idx - 1];
        if (idx % 2 == 0) {
            auto sk_left_j = party.left[j];   // shared_ptr<SecretKeyComponent>
            auto di = b.Didk_j(false, *sk_left_j, ID, j);
            shares.push_back(std::move(di));
        } else {
            auto sk_right_j = party.right[j]; // shared_ptr<SecretKeyComponent>
            auto di = b.Didk_j(true, *sk_right_j, ID, j);
            shares.push_back(std::move(di));
        }
    }

    auto idk = b.combIdk_j(J, shares, j);
    auto k_dec = b.dec(idk, ct);

    std::string k_enc_str = k_enc.getStr();
    std::string k_dec_str = k_dec.getStr();
    std::cout << "Encapsulated key: " << k_enc_str << "\n";
    std::cout << "Recovered   key: " << k_dec_str << "\n";

    bool equal = (k_enc == k_dec);
    std::cout << "Keys match? " << (equal ? "YES" : "NO") << "\n";
}

void testCCA_BT_KEM(int n, int t, int ell){
    using namespace ccakem;
    std::mt19937_64 rng{ std::random_device{}() };
    std::cout << "n=" << n << ", t=" << t << ", ell=" << ell << "\n";
    CCA_BT_KEM b = CCA_BT_KEM();

    auto keygenOut = b.KGen(n, t, ell, 128);

    std::uniform_int_distribution<int> distEll(0, ell - 1);
    int j = distEll(rng);

    auto encr = b.Encrypt(keygenOut.pk, j);

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

    std::vector<std::pair<bool, std::shared_ptr<btib::idkShare>>> shares;
    shares.reserve(J.size());
    for (int idx : J) {
        const auto& party = keygenOut.parties[idx - 1];
        if (idx % 2 == 0) {
            auto sk_left_j = party.left[j];   // shared_ptr<SecretKeyComponent>
            auto di = b.PDec(false, *sk_left_j, encr.second.first);
            shares.push_back(std::move(di));
        } else {
            auto sk_right_j = party.right[j]; // shared_ptr<SecretKeyComponent>
            auto di = b.PDec(true, *sk_right_j, encr.second.first);
            shares.push_back(std::move(di));
        }
    }

    auto k_dec = b.Comb(J, shares, encr.second.first);

    std::string k_enc_str = encr.first.getStr();
    std::string k_dec_str = k_dec.getStr();
    std::cout << "Encapsulated key: " << k_enc_str << "\n";
    std::cout << "Recovered   key: " << k_dec_str << "\n";

    bool equal = (encr.first == k_dec);
    std::cout << "Keys match? " << (equal ? "YES" : "NO") << "\n";
}

void testBTBF(int n, int t, int ell) {
    std::mt19937_64 rng{ std::random_device{}() };
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

void testTTT(int n, int t, int lambda, double e){
    using namespace ttt;
    std::mt19937_64 rng{ std::random_device{}() };
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

    auto D = [&J, &keygenOut, &E](const Ciphertext& c, const GT& k) {
        // Simulated decryption oracle
        std::vector<std::pair<bool, std::shared_ptr<idkShare>>> shares;
        shares.reserve(J.size());
        for (int idx : J) {
            shares.push_back(E.dec(keygenOut.parties[idx], idx, c));
        }
        auto k_rec = E.combine(J, shares, c);
        return k == k_rec;
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

int main(){
    using namespace ttt;
    using namespace btib;

    std::mt19937_64 rng{ std::random_device{}() };
    int ns[] = {100, 200};
    int ts[] = {25};

    int lambdas[] = {16};
    double es[] = {0.1};

    int ells[] = {10};

    for(int n : ns){
        for(int t : ts){

            for(int lambda : lambdas){
                for(double e : es){
                    //testTTT(n, t, lambda, e);
                }
            }

            for(int ell : ells){
                testCCA_BT_KEM(n, t, ell);
                //testBTIB(n, t, ell);
                //testBTBF(n, t, ell);
            }
        }
    }
    return 0;
}
