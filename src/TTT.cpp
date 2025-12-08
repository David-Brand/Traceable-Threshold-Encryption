#include "TTT.hpp"

namespace ttt{
    TTT::TTT(int t, int security_lambda) : t_(t), security_lambda_(security_lambda){
        fingerPrintingCode_ =  new fingerprinting::TardosCodes(t, 1 / std::pow(2, security_lambda));
        BTIB* tmp = new BTIB_KEM1();
        BTE = new CCA_BT_KEM(tmp);
    }

    TTTKeyGenOutput TTT::keygen(int n, int t, int security_lambda, double e){
        TTTKeyGenOutput out;
        float delta = (1/2 - e) / (1/2 - 2/std::sqrt(security_lambda));

        if(t != t_ || security_lambda != security_lambda_){
            t_ = t;
            security_lambda_ = security_lambda;
            fingerPrintingCode_ = new fingerprinting::TardosCodes(t, 1 / std::pow(2, security_lambda));
        }

        std::size_t l = fingerPrintingCode_->getLength();
        BTIBKeyGenOutput btOut = BTE->KGen(n, t, l, security_lambda);
        out.pk = btOut.pk;

        for(int i = 0; i<n; i++){
            SecretKey sk;
            auto [fingerprint, U] = fingerPrintingCode_->writeCodeWord();
            for(std::size_t j = 0; j<l; j++){
                if (fingerprint.get(j)){    //right key
                    sk.sk_components.push_back(btOut.parties[i].right[j]);
                } else{ //left key
                    sk.sk_components.push_back(btOut.parties[i].left[j]);
                }
            }
            out.parties.push_back(sk);
            out.tk.push_back(U);
        }
        return out;
    }

    std::pair<GT, std::pair<Ciphertext, proof>> TTT::enc(const std::shared_ptr<PublicKey>& pk){
        std::mt19937_64 random(std::random_device{}());
        std::uniform_int_distribution<int> dist(1, fingerPrintingCode_->getLength());
        return BTE->Encrypt(pk, dist(random));
    }

    std::pair<bool, std::shared_ptr<idkShare>> TTT::dec(const SecretKey& sk, int i, const Ciphertext& c){
        return BTE->PDec(fingerPrintingCode_->getBit(i, c.j), *sk.sk_components[c.j], c);
    }

    GT TTT::combine(const std::vector<int>& J, const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, const Ciphertext& c){
        return BTE->Comb(J, shares, c);
    }

    std::vector<std::size_t> TTT::trace(const std::shared_ptr<PublicKey>& pk, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext&, const GT&)>& D){
        auto N = std::pow(security_lambda_, 2);
        auto B = std::pow(security_lambda_, 3.0/2.0);

        auto totalLength = fingerPrintingCode_->getLength();

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
        return fingerprinting::TardosCodes::trace(tk, x, fingerPrintingCode_->Z());
    }

    int TTT::TrD(const std::shared_ptr<PublicKey>& pk, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext&, const GT&)>& D){
        int ctr = 0;

        #pragma omp parallel for reduction(+:ctr)
        for (int h = 0; h < (int)N; h++) {
            auto c0 = BTE->Encrypt(pk, j);
            auto c1 = BTE->Encrypt(pk, j);
            Ciphertext c;
            c.j = j;
            c.c0 = b0 ? c1.second.first.c0 : c0.second.first.c0;
            c.c1 = b1 ? c1.second.first.c1 : c0.second.first.c1;
            if(D(c, bk ? c1.first : c0.first)){
                ctr++;
            }
        }
        return ctr;
    }

}