// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: GPL-3.0-or-later

#include "btddh.hpp"
#include "Shamir.hpp"

using namespace mcl::bls12;

namespace btddh {
    bool BTDDH::initialized_ = false;
    G BTDDH::g_;

    // Initializes the BTDDH system for BLS12-381.
    // Sets up pairing, hash-to-curve mode, and defines the base group generator.
    void BTDDH::init(){
        // Prevent multiple initializations
        if(initialized_) return;
        // Initialize pairing and hash-to-curve
        initPairing(mcl::BLS12_381);
        setMapToMode(MCL_MAP_TO_MODE_HASH_TO_CURVE);
        // Load fixed generator for G
        const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
        g_.setStr(g1Str, 16);

        initialized_ = true;
    }

    KeyGenOutput_ddh BTDDH::keygen(int num_parties, int decryption_threshold, int fingerprint_length, int security_lambda){
        KeyGenOutput_ddh output;
        output.public_key = new PublicKey_ddh();
        output.n = num_parties;
        output.t = decryption_threshold;
        output.ell = fingerprint_length;
        output.secret_keys.resize(num_parties);

        for(int j=0; j<fingerprint_length; j++){
            //sample random x,y,z in Zp
            Zp x, y, z;
            x.setByCSPRNG();
            y.setByCSPRNG();
            z.setByCSPRNG();

            // calculate public key elements and put them into output
            G X, Y, Z;
            G::mul(X, g_, x);
            G::mul(Y, g_, y);
            G::mul(Z, g_, z);
            output.public_key->X.push_back(X);
            output.public_key->Y.push_back(Y);
            output.public_key->Z.push_back(Z);

            // generate secret key components for each party
            std::vector<Zp> s;
            s = share<Zp>(x, num_parties, decryption_threshold);
            for(int i=0; i<num_parties; i++){
                SecretKeyComponent_ddh left, right;
                Zp::div(left.secret, s[i], y); // sk_i,0^{(j)} = s_i/y
                Zp::div(right.secret, s[i], z); // sk_i,1^{(j)} = s_i/z
                output.secret_keys[i].right.push_back(right);
                output.secret_keys[i].left.push_back(left);
            }
        }
        return output;
    }

    std::pair<G, Ciphertext_ddh> BTDDH::encaps(const PublicKey_ddh* public_key, int j){
        // sample random r in Zp
        Zp r;
        r.setByCSPRNG();

        // compute k = X_j^{r}
        G k;
        G::mul(k, public_key->X[j], r);

        // compute ciphertext components
        Ciphertext_ddh c;
        c.j = j;
        G::mul(c.c0.comp, public_key->Y[j-1], r); // Y_j^r
        G::mul(c.c1.comp, public_key->Z[j-1], r); // Z_j^r

        return {k, c};
    }

    G BTDDH::decShare(const SecretKeyComponent_ddh& sk_i_b, const CiphertextComponent_ddh& c_b){
        G d_i;
        G::mul(d_i, c_b.comp, sk_i_b.secret); // d_i = c_b^{sk_i,b^{(j)}}
        return d_i;
    }

    G BTDDH::combine(const std::vector<int>& providing_parties, const std::vector<G>& shares){
        G k;

        for(int i=0; i<providing_parties.size(); i++){
            int idx_i = providing_parties[i];
            // compute Lagrange coefficient
            Zp lambda = lagrangeCoeff(providing_parties, idx_i);
            // k += d_i^{lambda}
            G temp;
            G::mul(temp, shares[i], lambda);
            G::add(k, k, temp);
        }
        return k;
    }
}
