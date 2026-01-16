#pragma once

// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: MIT

#include <mcl/bls12_381.hpp>
#include <mcl/curve_type.h>

#include <array>
#include <cstdint>
#include <utility>
#include <vector>
#include <string>

namespace btbf {
    using G1 = mcl::bls12::G1;
    using G2 = mcl::bls12::G2;
    using GT = mcl::bls12::GT;
    using Fr = mcl::bls12::Fr;

    // H2 key
    using SymKey = std::vector<uint8_t>;

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
        G1 Q;
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
        PublicKey* public_key;
        // pkc = ⊥ in the scheme --> omitted
        std::vector<PartySecret> parties; // size num_parties
        int n;
        int t;
        int ell;
    };

class BTBF {
public:
    // Must be called once before using any other BTBF functions.
    static void init();

    // BTBF.KeyGen(1^λ, num_parties, decryption_threshold, ℓ)
    static KeyGenOutput keygen(int num_parties, int decryption_threshold, int fingerprint_length, int security_lambda);

    // BTBF.Enc(public_key, j): returns (k, c)
    static std::pair<SymKey, Ciphertext> encaps(const PublicKey* public_key, int j);

    // BTBF.Dec(j, sk_i,b^{(j)}, c_b) -> d_i in GT
    static GT decShare(const SecretKeyComponent& sk_i_b, const CiphertextComponent& c_b);

    // BTBF.Combine(pkc = ⊥, j, c, providing_parties, {d_i}_{i∈providing_parties}) -> k
    // providing_parties : 1-based indices of parties; shares : d_i in the same order as providing_parties.
    static SymKey combine(const std::vector<int>& providing_parties, const std::vector<GT>& shares);

    static int getLambda(){
        return security_lambda_bits_;
    }

    // H1 : N -> G1
    static G1 H1(int j);

private:
    static bool initialized_;
    static G1 g1_; // generator of G1
    static G2 g2_; // generator of G2
    static int security_lambda_bits_;

    // H2 : GT -> {0,1}^λ
    static SymKey H2(const GT& w);

    //Lagrange coefficient λ_i^providing_parties for point i in providing_parties
    static Fr lagrangeCoeff(const std::vector<int>& J, int i);
};

} // namespace btbf
