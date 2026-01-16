// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: MIT

#include "TTT.hpp"

namespace ttt{
    // Constructs the BTBF-based TTT wrapper with given threshold/security.
    // Initializes the fingerprinting scheme and the BTBF backend.
    TTT_BTBF::TTT_BTBF(int decryption_threshold, int security_lambda) : decryption_threshold_(decryption_threshold), security_lambda_(security_lambda){
        // Initialize tracing code parameters
        fingerPrintingCode_ =  new fingerprinting::TardosCodes(decryption_threshold, 1 / std::pow(2, security_lambda));
        // Initialize encryption backend
        BTE = new BTBF();
        BTE->init();
    }

    // Generates keys for all parties using BTBF and the fingerprinting scheme.
    // Selects left/right secret key components per user based on their fingerprint.
    TTTKeyGenOutput TTT_BTBF::keygen(int num_parties, int decryption_threshold, int security_lambda, double decoder_error){
        TTTKeyGenOutput out;
        float delta = (1/2 - decoder_error) / (1/2 - 2/std::sqrt(security_lambda));

        // Refresh internal state if parameters changed
        if(decryption_threshold != decryption_threshold_ || security_lambda != security_lambda_){
            decryption_threshold_ = decryption_threshold;
            security_lambda_ = security_lambda;
            fingerPrintingCode_ = new fingerprinting::TardosCodes(decryption_threshold, 0.1);
        }

        // Generate BTBF keys and map them via fingerprints
        std::size_t fingerprint_length = fingerPrintingCode_->getLength();
        KeyGenOutput btOut = BTE->keygen(num_parties, decryption_threshold, fingerprint_length, security_lambda);
        out.public_key = btOut.public_key;

        // Assign per-user secret components according to fingerprint bits
        for(int i = 0; i<num_parties; i++){
            SecretKey sk;
            auto [fingerprint, U] = fingerPrintingCode_->writeCodeWord();
            for(std::size_t j = 0; j<fingerprint_length; j++){
                if (fingerprint.get(j)){    //right key
                    sk.sk_components.push_back(btOut.parties[i].right[j]);
                } else{ //left key
                    sk.sk_components.push_back(btOut.parties[i].left[j]);
                }
            }
            // Store user keys and tracing info
            out.parties.push_back(sk);
            out.tracing_key.push_back(U);
        }
        return out;
    }

    // Encapsulates to a random fingerprint index.
    // Delegates to BTBF encapsulation.
    std::pair<SymKey, Ciphertext> TTT_BTBF::enc(const PublicKey* public_key){
        // Pick random fingerprint position
        std::mt19937_64 random(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, fingerPrintingCode_->getLength());
        return BTE->encaps(public_key, dist(random));
    }

    // Produces a decryption share for party i based on their selected component.
    // Chooses c0 or c1 depending on the fingerprint bit.
    GT TTT_BTBF::dec(const SecretKey& sk, int i, const Ciphertext& c){
        // Select component by fingerprint bit
        if(fingerPrintingCode_->getBit(i, c.j)){
            return BTE->decShare(sk.sk_components[c.j], c.c0);
        } else {
            return BTE->decShare(sk.sk_components[c.j], c.c1);
        }
    }

    // Combines decryption shares using BTBF to recover the symmetric key.
    SymKey TTT_BTBF::combine(const std::vector<int>& J, const std::vector<GT>& shares){
        // Delegate combination to backend
        return BTE->combine(J, shares);
    }

    // Traces colluders by reconstructing a bitset via decoder oracles and Tardos tracing.
    // Uses Monte Carlo sampling (TrD) to decide bit values and then accuses via U-weights.
    std::vector<std::size_t> TTT_BTBF::trace(const PublicKey* public_key, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext&, const SymKey&)>& D){
        // Choose sampling and threshold parameters
        auto N = std::pow(security_lambda_, 2);
        auto B = std::pow(security_lambda_, 3.0/2.0);

        auto totalLength = fingerPrintingCode_->getLength();

        // Recover pirate word bits by oracle comparisons
        PackedBitset x(totalLength);

        for(std::size_t j=1; j<totalLength+1; j++){
            auto p001 = TrD(public_key, j, N, 0, 0, 1, D);
            auto p100 = TrD(public_key, j, N, 1, 0, 0, D);
            auto p111 = TrD(public_key, j, N, 1, 1, 1, D);

            auto a0 = std::abs(p001 - p100);
            auto a1 = std::abs(p001 - p111);
            if(a0 >= B){
                x.set(j, false);
            } else if (a1 >= B){
                x.set(j, true);
            } else {
                x.set(j, false);
            }
        }
        // Accuse via Tardos tracing
        return fingerprinting::TardosCodes::trace(tk, x, fingerPrintingCode_->get_accusation_threshold());
    }

    // Runs the decoder oracle N times with selected bit choices and counts accepts.
    // Builds mixed ciphertexts by toggling components and compares with key choices.
    int TTT_BTBF::TrD(const PublicKey* public_key, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext&, const SymKey&)>& D){
        int ctr = 0;

        // Monte Carlo sampling loop
        #pragma omp parallel for reduction(+:ctr)
        for (int h = 0; h < (int)N; h++) {
            // Create two encapsulations and mix components
            auto c0 = BTE->encaps(public_key, j);
            auto c1 = BTE->encaps(public_key, j);
            Ciphertext c;
            c.j = j;
            c.c0 = b0 ? c1.second.c0 : c0.second.c0;
            c.c1 = b1 ? c1.second.c1 : c0.second.c1;
            // Query decoder oracle and count success
            if(D(c, bk ? c1.first : c0.first)){
                ctr++;
            }
        }
        return ctr;
    }

//-----------------------------------------------TTT_BTDDH----------------------------------------------------

    // Constructs the BTDDH-based TTT wrapper and initializes backend.
    TTT_BTDDH::TTT_BTDDH(int decryption_threshold, int security_lambda) : decryption_threshold_(decryption_threshold), security_lambda_(security_lambda){
        // Initialize tracing code and backend
        fingerPrintingCode_ =  new fingerprinting::TardosCodes(decryption_threshold, 1 / std::pow(2, security_lambda));
        BTE = new BTDDH();
        BTE->init();
    }

    // Generates BTDDH keys mapped via fingerprints to left/right components.
    TTTKeyGenOutput_ddh TTT_BTDDH::keygen(int num_parties, int decryption_threshold, int security_lambda, double decoder_error){
        TTTKeyGenOutput_ddh out;
        float delta = (1/2 - decoder_error) / (1/2 - 2/std::sqrt(security_lambda));

        // Refresh code if parameters change
        if(decryption_threshold != decryption_threshold_ || security_lambda != security_lambda_){
            decryption_threshold_ = decryption_threshold;
            security_lambda_ = security_lambda;
            fingerPrintingCode_ = new fingerprinting::TardosCodes(decryption_threshold, 0.1);
        }

        // Generate base keys and assemble per-user secrets
        std::size_t fingerprint_length = fingerPrintingCode_->getLength();
        KeyGenOutput_ddh btOut = BTE->keygen(num_parties, decryption_threshold, fingerprint_length, security_lambda);
        out.public_key = btOut.public_key;

        // Assign left/right secrets per fingerprint bit
        for(int i = 0; i<num_parties; i++){
            SecretKey_ddh sk;
            auto [fingerprint, U] = fingerPrintingCode_->writeCodeWord();
            for(std::size_t j = 0; j<fingerprint_length; j++){
                if (fingerprint.get(j)){    //right key
                    sk.sk_components.push_back(btOut.secret_keys[i].right[j]);
                } else{ //left key
                    sk.sk_components.push_back(btOut.secret_keys[i].left[j]);
                }
            }
            out.parties.push_back(sk);
            out.tracing_key.push_back(U);
        }
        return out;
    }

    // Encapsulates to a random fingerprint index using BTDDH.
    std::pair<G, Ciphertext_ddh> TTT_BTDDH::enc(const PublicKey_ddh* public_key){
        // Pick random fingerprint position
        std::mt19937_64 random(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, fingerPrintingCode_->getLength());
        return BTE->encaps(public_key, dist(random));
    }

    // Produces a decryption share for party i based on the chosen component.
    G TTT_BTDDH::dec(const SecretKey_ddh& sk, int i, const Ciphertext_ddh& c){
        // Select component by fingerprint bit
        if(fingerPrintingCode_->getBit(i, c.j)){
            return BTE->decShare(sk.sk_components[c.j], c.c0);
        } else {
            return BTE->decShare(sk.sk_components[c.j], c.c1);
        }
    }

    // Combines shares to recover the group element key in BTDDH.
    G TTT_BTDDH::combine(const std::vector<int>& J, const std::vector<G>& shares){
        // Delegate to backend combination
        return BTE->combine(J, shares);
    }

    // Traces BTDDH colluders by reconstructing pirate bits via TrD and Tardos tracing.
    std::vector<std::size_t> TTT_BTDDH::trace(const PublicKey_ddh* public_key, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext_ddh&, const G&)>& D){
        // Choose sampling and threshold parameters
        auto N = std::pow(security_lambda_, 2);
        auto B = std::pow(security_lambda_, 3.0/2.0);

        auto totalLength = fingerPrintingCode_->getLength();

        // Recover pirate word bits by oracle comparisons
        PackedBitset x(totalLength);

        for(std::size_t j=1; j<totalLength+1; j++){
            auto p001 = TrD(public_key, j, N, 0, 0, 1, D);
            auto p100 = TrD(public_key, j, N, 1, 0, 0, D);
            auto p111 = TrD(public_key, j, N, 1, 1, 1, D);

            auto a0 = std::abs(p001 - p100);
            auto a1 = std::abs(p001 - p111);
            if(a0 >= B){
                x.set(j, false);
            } else if (a1 >= B){
                x.set(j, true);
            } else {
                x.set(j, false);
            }
        }
        // Accuse via Tardos tracing
        return fingerprinting::TardosCodes::trace(tk, x, fingerPrintingCode_->get_accusation_threshold());
    }

    // Runs the BTDDH decoder oracle N times and counts accepts for bit decisions.
    int TTT_BTDDH::TrD(const PublicKey_ddh* public_key, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext_ddh&, const G&)>& D){
        int ctr = 0;

        // Monte Carlo sampling loop
        #pragma omp parallel for reduction(+:ctr)
        for (int h = 0; h < (int)N; h++) {
            // Create two encapsulations and mix components
            auto c0 = BTE->encaps(public_key, j);
            auto c1 = BTE->encaps(public_key, j);
            Ciphertext_ddh c;
            c.j = j;
            c.c0 = b0 ? c1.second.c0 : c0.second.c0;
            c.c1 = b1 ? c1.second.c1 : c0.second.c1;
            // Query decoder oracle and count success
            if(D(c, bk ? c1.first : c0.first)){
                ctr++;
            }
        }
        return ctr;
    }
}