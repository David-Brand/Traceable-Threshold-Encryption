#include "TTT.hpp"

namespace ttt{
    TTT_BTBF::TTT_BTBF(int decryption_threshold, int security_lambda) : decryption_threshold_(decryption_threshold), security_lambda_(security_lambda){
        fingerPrintingCode_ =  new fingerprinting::TardosCodes(decryption_threshold, 1 / std::pow(2, security_lambda));
        BTE = new BTBF();
        BTE->init();
    }

    TTTKeyGenOutput TTT_BTBF::keygen(int num_parties, int decryption_threshold, int security_lambda, double decoder_error){
        TTTKeyGenOutput out;
        float delta = (1/2 - decoder_error) / (1/2 - 2/std::sqrt(security_lambda));

        if(decryption_threshold != decryption_threshold_ || security_lambda != security_lambda_){
            decryption_threshold_ = decryption_threshold;
            security_lambda_ = security_lambda;
            fingerPrintingCode_ = new fingerprinting::TardosCodes(decryption_threshold, 0.1);
        }

        std::size_t fingerprint_length = fingerPrintingCode_->getLength();
        KeyGenOutput btOut = BTE->keygen(num_parties, decryption_threshold, fingerprint_length, security_lambda);
        out.public_key = btOut.public_key;

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
            out.parties.push_back(sk);
            out.tracing_key.push_back(U);
        }
        return out;
    }

    std::pair<SymKey, Ciphertext> TTT_BTBF::enc(const PublicKey* public_key){
        std::mt19937_64 random(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, fingerPrintingCode_->getLength());
        return BTE->encaps(public_key, dist(random));
    }

    GT TTT_BTBF::dec(const SecretKey& sk, int i, const Ciphertext& c){
        if(fingerPrintingCode_->getBit(i, c.j)){
            return BTE->decShare(sk.sk_components[c.j], c.c0);
        } else {
            return BTE->decShare(sk.sk_components[c.j], c.c1);
        }
    }

    SymKey TTT_BTBF::combine(const std::vector<int>& J, const std::vector<GT>& shares){
        return BTE->combine(J, shares);
    }

    std::vector<std::size_t> TTT_BTBF::trace(const PublicKey* public_key, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext&, const SymKey&)>& D){
        auto N = std::pow(security_lambda_, 2);
        auto B = std::pow(security_lambda_, 3.0/2.0);

        auto totalLength = fingerPrintingCode_->getLength();

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
        return fingerprinting::TardosCodes::trace(tk, x, fingerPrintingCode_->get_accusation_threshold());
    }

    int TTT_BTBF::TrD(const PublicKey* public_key, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext&, const SymKey&)>& D){
        int ctr = 0;

        #pragma omp parallel for reduction(+:ctr)
        for (int h = 0; h < (int)N; h++) {
            auto c0 = BTE->encaps(public_key, j);
            auto c1 = BTE->encaps(public_key, j);
            Ciphertext c;
            c.j = j;
            c.c0 = b0 ? c1.second.c0 : c0.second.c0;
            c.c1 = b1 ? c1.second.c1 : c0.second.c1;
            if(D(c, bk ? c1.first : c0.first)){
                ctr++;
            }
        }
        return ctr;
    }

//-----------------------------------------------TTT_BTDDH----------------------------------------------------

    TTT_BTDDH::TTT_BTDDH(int decryption_threshold, int security_lambda) : decryption_threshold_(decryption_threshold), security_lambda_(security_lambda){
        fingerPrintingCode_ =  new fingerprinting::TardosCodes(decryption_threshold, 1 / std::pow(2, security_lambda));
        BTE = new BTDDH();
        BTE->init();
    }

    TTTKeyGenOutput_ddh TTT_BTDDH::keygen(int num_parties, int decryption_threshold, int security_lambda, double decoder_error){
        TTTKeyGenOutput_ddh out;
        float delta = (1/2 - decoder_error) / (1/2 - 2/std::sqrt(security_lambda));

        if(decryption_threshold != decryption_threshold_ || security_lambda != security_lambda_){
            decryption_threshold_ = decryption_threshold;
            security_lambda_ = security_lambda;
            fingerPrintingCode_ = new fingerprinting::TardosCodes(decryption_threshold, 0.1);
        }

        std::size_t fingerprint_length = fingerPrintingCode_->getLength();
        KeyGenOutput_ddh btOut = BTE->keygen(num_parties, decryption_threshold, fingerprint_length, security_lambda);
        out.public_key = btOut.public_key;

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

    std::pair<G, Ciphertext_ddh> TTT_BTDDH::enc(const PublicKey_ddh* public_key){
        std::mt19937_64 random(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, fingerPrintingCode_->getLength());
        return BTE->encaps(public_key, dist(random));
    }

    G TTT_BTDDH::dec(const SecretKey_ddh& sk, int i, const Ciphertext_ddh& c){
        if(fingerPrintingCode_->getBit(i, c.j)){
            return BTE->decShare(sk.sk_components[c.j], c.c0);
        } else {
            return BTE->decShare(sk.sk_components[c.j], c.c1);
        }
    }

    G TTT_BTDDH::combine(const std::vector<int>& J, const std::vector<G>& shares){
        return BTE->combine(J, shares);
    }

    std::vector<std::size_t> TTT_BTDDH::trace(const PublicKey_ddh* public_key, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext_ddh&, const G&)>& D){
        auto N = std::pow(security_lambda_, 2);
        auto B = std::pow(security_lambda_, 3.0/2.0);

        auto totalLength = fingerPrintingCode_->getLength();

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
        return fingerprinting::TardosCodes::trace(tk, x, fingerPrintingCode_->get_accusation_threshold());
    }

    int TTT_BTDDH::TrD(const PublicKey_ddh* public_key, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext_ddh&, const G&)>& D){
        int ctr = 0;

        #pragma omp parallel for reduction(+:ctr)
        for (int h = 0; h < (int)N; h++) {
            auto c0 = BTE->encaps(public_key, j);
            auto c1 = BTE->encaps(public_key, j);
            Ciphertext_ddh c;
            c.j = j;
            c.c0 = b0 ? c1.second.c0 : c0.second.c0;
            c.c1 = b1 ? c1.second.c1 : c0.second.c1;
            if(D(c, bk ? c1.first : c0.first)){
                ctr++;
            }
        }
        return ctr;
    }
}