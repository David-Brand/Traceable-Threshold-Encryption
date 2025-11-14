#pragma once
#include <mcl/bls12_381.hpp>
#include <string>
#include <vector>
#include <utility>

namespace btbf {

using Fr = mcl::bls12::Fr;
using G1 = mcl::bls12::G1;
using G2 = mcl::bls12::G2;
using GT = mcl::bls12::GT;

// Public key: X = g1^(αyz), Y = g1^y, Z = g1^z
struct PublicKey {
    G1 X, Y, Z;
};

// A single (k0, k1) secret share: k0 ∈ G1, k1 ∈ G2
struct Share {
    G1 k0;
    G2 k1;
};

using LRShare = std::pair<Share, Share>;

// Secret key for one party: ℓ positions of (left, right) shares
struct PartySecretKey {
    std::vector<LRShare> pos; // size ℓ
};

// Ciphertext side c_b = (u_b ∈ G2, v_b ∈ G1)
struct Side {
    G2 u;
    G1 v;
};

struct Ciphertext {
    Side c0;
    Side c1;
};

//output key
struct SharedKey {
    std::string hex;
};

class BTBF {
public:
    //curve setup
    static void init();

    // Key generation
    static void KeyGen(int n, int t, int ell, PublicKey& pk, std::vector<PartySecretKey>& sks);

    // Encapsulation to position j (1-based): returns (c, k)
    static void Enc(const PublicKey& pk, int j, Ciphertext& ct, SharedKey& k);

    // A party's decryption share using its left/right bit b for position j (1-based)
    static GT DecShare(const PartySecretKey& ski, int j, int b, const Ciphertext& ct);

    static SharedKey Combine(const std::vector<int>& J, const std::vector<GT>& shares);

private:
    //N -> G1
    static void H1(G1& out, int j);

    static void G1gen(G1& g1);
    static void G2gen(G2& g2);

    //GT -> {0,1}^λ
    static std::string kdf(const GT& x);

    // Lagrange coefficients at 0 for indices J (in Fr)
    static std::vector<Fr> lagrangeAtZero(const std::vector<int>& J);
};

} // namespace btbf
