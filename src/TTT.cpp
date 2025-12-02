#include "TTT.hpp"

namespace ttt{
    auto tmp = fingerprinting::TardosCodes(2, 0.5);

    TTT::TTT(int t, int security_lambda) : t_(t), security_lambda_(security_lambda), fingerPrintingCode_(tmp) {
        fingerPrintingCode_ = fingerprinting::TardosCodes(t, 1 / std::pow(2, security_lambda));
        BTBF::init();
    }

    TTTKeyGenOutput TTT::keygen(int n, int t, int security_lambda, double e) const{
        TTTKeyGenOutput out;
        float delta = (1/2 - e) / (1/2 - 2/std::sqrt(security_lambda));

        if(t != t_ || security_lambda != security_lambda_){
            fingerPrintingCode_ = fingerprinting::TardosCodes(t, 1 / std::pow(2, security_lambda));
        }

        std::size_t l = fingerPrintingCode_.getLength();

        BTBF::KeyGenOutput btbfOut = BTBF::keygen(n, t, l, security_lambda);
        out.pk = btbfOut.pk;

        for(int i = 0; i<n; i++){
            std::vector<SecretKeyComponent_b> sk_b_vec;
            auto [fingerprint, U] = fingerprinting::TardosCodes::generateCodeWord(fingerPrintingCode_.getProbabilities());
            for(std::size_t j = 0; j<l; j++){
                SecretKeyComponent_b sk;
                if (fingerprint.get(j)){    //right key
                    sk.sk_b = btbfOut.parties[i].right[j];
                    sk.b = true;
                } else{ //left key
                    sk.sk_b = btbfOut.parties[i].left[j];
                    sk.b = false;
                }
                sk_b_vec.push_back(sk);
            }
            out.parties.push_back(sk_b_vec);
            out.tk.push_back(U);
        }
        return out;
    }

    std::pair<BTBF::SymKey, BTBF::Ciphertext> TTT::enc(const BTBF::PublicKey& pk) const{
        std::mt19937_64 random(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, fingerPrintingCode_.getLength());
        return BTBF::encaps(pk, dist(random));
    }

    TTT::GT TTT::dec(const SecretKeyComponent_b& sk_b, const BTBF::Ciphertext& c){
        return BTBF::decShare(sk_b.sk_b, sk_b.b ? c.c1 : c.c0);
    }

    BTBF::SymKey TTT::combine(const std::vector<int>& J, const std::vector<GT>& shares){
        return BTBF::combine(J, shares);
    }

    std::vector<std::size_t> TTT::trace(BTBF::PublicKey pk, std::vector<std::vector<double>> tk, const std::function<bool(const BTBF::Ciphertext&, const BTBF::SymKey&)>& D) const{
        int lambda = BTBF::getLambda();
        auto N = std::pow(lambda, 2);
        auto B = std::pow(lambda, 3.0/2.0);

        auto totalLength = fingerPrintingCode_.getLength();

        PackedBitset x(totalLength);

        for(std::size_t j=1; j<totalLength+1; j++){
            auto p001 = TrD(pk, j, N, 0, 0, 1, D);
            auto p100 = TrD(pk, j, N, 1, 0, 0, D);
            auto p111 = TrD(pk, j, N, 1, 1, 1, D);

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
        return fingerprinting::TardosCodes::trace(tk, x, fingerPrintingCode_.Z());
    }

    int TTT::TrD(BTBF::PublicKey pk, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const BTBF::Ciphertext&, const BTBF::SymKey&)>& D) const{
        int ctr = 0;

        #pragma omp parallel for reduction(+:ctr)
        for (int h = 0; h < (int)N; h++) {
            auto c0 = BTBF::encaps(pk, j);
            auto c1 = BTBF::encaps(pk, j);
            BTBF::Ciphertext c;
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