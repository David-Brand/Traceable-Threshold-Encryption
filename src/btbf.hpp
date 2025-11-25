#pragma once

#include <mcl/bls12_381.hpp>
#include <mcl/curve_type.h>

#include <array>
#include <cstdint>
#include <utility>
#include <vector>
#include <string>

namespace btbf {

class BTBF {
public:
    using G1 = mcl::bn::G1;
    using G2 = mcl::bn::G2;
    using GT = mcl::bn::GT;
    using Fr = mcl::bn::Fr;

    // H2 key = 256-bit string
    using SymKey = std::array<uint8_t, 32>;

    struct SecretKeyComponent {
        G1 k0; // in G1
        G2 k1; // in G2
    };

    // For each party i, for each position j we have left (b=0) and right (b=1) keys
    struct PartySecret {
        std::vector<SecretKeyComponent> left;   // sk_i,0^{(j)} for j=1..ℓ
        std::vector<SecretKeyComponent> right;  // sk_i,1^{(j)} for j=1..ℓ
    };

    struct PublicKey {
        G1 X; // g1^{α y z}
        G1 Y; // g1^y
        G1 Z; // g1^z
    };

    struct CiphertextComponent {
        G2 u; // in G2
        G1 v; // in G1
    };

    struct Ciphertext {
        CiphertextComponent c0; // encrypts Y^r
        CiphertextComponent c1; // encrypts Z^r
        int j;                  // position j this ciphertext is for
    };

    struct KeyGenOutput {
        PublicKey pk;
        // pkc = ⊥ in the scheme --> omitted
        std::vector<PartySecret> parties; // size n
        int n;
        int t;
        int ell;
    };

    // --- API ---

    // Must be called once before using any other BTBF functions.
    static void init();

    // BTBF.KeyGen(1^λ, n, t, ℓ)
    static KeyGenOutput keygen(int n, int t, int ell);

    // BTBF.Enc(pk, j): returns (k, c)
    static std::pair<SymKey, Ciphertext> encaps(const PublicKey& pk, int j);

    // BTBF.Dec(j, sk_i,b^{(j)}, c_b) -> d_i in GT
    static GT decShare(const SecretKeyComponent& sk_i_b,
                       const CiphertextComponent& c_b);

    // BTBF.Combine(pkc = ⊥, j, c, J, {d_i}_{i∈J}) -> k
    // J : 1-based indices of parties; shares : d_i in the same order as J.
    static SymKey combine(const std::vector<int>& J,
                          const std::vector<GT>& shares);

private:
    static bool initialized_;
    static G1 g1_; // "generator" used for X,Y,Z
    static G2 g2_; // "generator" used wherever g2 appears

    // H1 : N -> G1
    static G1 H1(int j);

    // H2 : GT -> {0,1}^λ  (here λ = 256)
    static SymKey H2(const GT& w);

    // helper : Lagrange coefficient λ_i^J at 0 for point i in J
    static Fr lagrangeCoeffAtZero(const std::vector<int>& J, int i);
};

} // namespace btbf
