#pragma once

// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: GPL-3.0-or-later

#include <mcl/bls12_381.hpp>
#include <mcl/curve_type.h>

#include <array>
#include <cstdint>
#include <utility>
#include <vector>
#include <string>

namespace btddh{
    using G = mcl::bls12::G1;
    using Zp = mcl::bls12::Fr;

    struct SecretKeyComponent_ddh{
        Zp secret;
    };

    struct PartySecret_ddh {
        std::vector<SecretKeyComponent_ddh> left;   // sk_i,0^{(j)} for j=1..ℓ
        std::vector<SecretKeyComponent_ddh> right;  // sk_i,1^{(j)} for j=1..ℓ
    };

    struct PublicKey_ddh {
        G Q;
        std::vector<G> X; // g1^{α y z}
        std::vector<G> Y; // g1^y
        std::vector<G> Z; // g1^z
    };

    struct CiphertextComponent_ddh {
        G comp; // in G
    };

    struct Ciphertext_ddh {
        CiphertextComponent_ddh c0; // encrypts Y^r
        CiphertextComponent_ddh c1; // encrypts Z^r
        int j;                  // position j this ciphertext is for
    };

    struct KeyGenOutput_ddh {
        PublicKey_ddh* public_key;
        // pkc = ⊥ in the scheme --> omitted
        std::vector<PartySecret_ddh> secret_keys; // size num_parties
        int n;
        int t;
        int ell;
    };

    class BTDDH{
    public:
        // Must be called once before using any other BTDDH functions.
        static void init();

        // BTDDH.KeyGen(1^λ, num_parties, decryption_threshold, ℓ)
        static KeyGenOutput_ddh keygen(int num_parties, int decryption_threshold, int fingerprint_length, int security_lambda);

        // BTDDH.Enc(public_key, j): returns (k, c)
        static std::pair<G, Ciphertext_ddh> encaps(const PublicKey_ddh* public_key, int j);

        // BTDDH.Dec(j, sk_i,b^{(j)}, c_b) -> d_i in G
        static G decShare(const SecretKeyComponent_ddh& sk_i_b, const CiphertextComponent_ddh& c_b);

        // BTDDH.Combine(pkc = ⊥, j, c, providing_parties, {d_i}_{i∈providing_parties}) -> k
        // providing_parties : 1-based indices of parties; shares : d_i in the same order as providing_parties.
        static G combine(const std::vector<int>& providing_parties, const std::vector<G>& shares);

    private:
        static bool initialized_;
        static G g_;

        static Zp lagrangeCoeff(const std::vector<int>& J, int i){
            Zp fi = i;
            Zp lambda = 1;

            for (int v : J) {
                if (v == i) continue;

                Zp fv = v;
                Zp sub;
                Zp::sub(sub, fv, fi);
                Zp div;
                Zp::div(div, fv, sub);
                Zp::mul(lambda, lambda, div);
            }
            return lambda;
        }
    };
}
