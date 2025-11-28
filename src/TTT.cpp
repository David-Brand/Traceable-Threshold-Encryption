#include "TTT.hpp"

namespace TTT{
    fingerprinting::LogLengthCodes defaultFingerPrintingCode(2, 1, 0.7, 0.1);

    TTT::TTT() : fingerPrintingCode_(defaultFingerPrintingCode) {
        BTBF::init();
    }

    TTTKeyGenOutput TTT::keygen(int n, int t, int security_lambda, double e) const{
        TTTKeyGenOutput out;
        float delta = (1/2 - e) / (1/2 - 2/std::sqrt(security_lambda));

        fingerPrintingCode_ = fingerprinting::LogLengthCodes(n, t, 1 / std::pow(2, security_lambda), delta);
        std::size_t l = fingerPrintingCode_.totalLength();
        std::size_t blockLength = fingerPrintingCode_.blockLength();

        BTBF::KeyGenOutput btbfOut = BTBF::keygen(n, t, l, security_lambda);
        out.pk = btbfOut.pk;

        for(int i = 0; i<n; i++){
            std::vector<SecretKeyComponent_b> sk_b_vec;
            std::vector<PackedBitset> fingerprint;
            fingerPrintingCode_.getLLCodeword(i, fingerprint);
            for(std::size_t j = 0; j<l; j++){
                SecretKeyComponent_b sk;
                if (fingerprint[std::floor(j / blockLength)].get(j % blockLength)){    //right key
                    sk.sk_b = btbfOut.parties[i].right[j];
                    sk.b = true;
                } else{ //left key
                    sk.sk_b = btbfOut.parties[i].left[j];
                    sk.b = false;
                }
                sk_b_vec.push_back(sk);
            }
            out.parties.push_back(sk_b_vec);
        }
        return out;
    }

    std::pair<BTBF::SymKey, BTBF::Ciphertext> TTT::enc(const BTBF::PublicKey& pk) const{
        return BTBF::encaps(pk, std::rand() % fingerPrintingCode_.totalLength());
    }

    TTT::GT TTT::dec(const SecretKeyComponent_b& sk_b, const BTBF::Ciphertext& c) const{
        return BTBF::decShare(sk_b.sk_b, sk_b.b ? c.c1 : c.c0);
    }

    BTBF::SymKey TTT::combine(const std::vector<int>& J, const std::vector<GT>& shares) const{
        return BTBF::combine(J, shares);
    }

    std::size_t TTT::trace(BTBF::PublicKey pk, const std::function<bool(const BTBF::Ciphertext&, const BTBF::SymKey&)>& D) const{
        int lambda = BTBF::getLambda();
        auto N = std::pow(lambda, 2);
        auto B = std::pow(lambda, 3.0/2.0);

        auto totalLength = fingerPrintingCode_.totalLength();
        auto blockLength = fingerPrintingCode_.blockLength();

        std::vector<PackedBitset> x(totalLength, PackedBitset(blockLength));
        std::vector<PackedBitset> x_unreadable(totalLength, PackedBitset(blockLength));

        for(std::size_t j=0; j<totalLength; j++){
            auto p001 = TrD(pk, j, N, 0, 0, 1, D);
            auto p100 = TrD(pk, j, N, 1, 0, 0, D);
            auto p111 = TrD(pk, j, N, 1, 1, 1, D);

            auto a0 = std::abs(p001 - p100);
            auto a1 = std::abs(p001 - p111);
            if(a0 >= B){
                x[std::floor(j / blockLength)].set(j%blockLength, false);
                x_unreadable[std::floor(j / blockLength)].set(j%blockLength, false);
            } else if (a1 >= B){
                x[std::floor(j / blockLength)].set(j%blockLength, true);
                x_unreadable[std::floor(j / blockLength)].set(j%blockLength, false);
            } else {
                x_unreadable[std::floor(j / blockLength)].set(j%blockLength, true);
                x[std::floor(j / blockLength)].set(j%blockLength, false);
            }
        }
        return fingerPrintingCode_.trace(x, x_unreadable);
    }

    int TTT::TrD(BTBF::PublicKey pk, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const BTBF::Ciphertext&, const BTBF::SymKey&)>& D) const{
        int ctr = 0;

        for(uint32_t i=0; i<N; i++){
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