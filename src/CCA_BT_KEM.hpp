#pragma once

#include "BTIB.hpp"
#include "SignatureScheme.hpp"
#include <memory>

namespace ccakem{
    using namespace btib;
    using namespace proofsystem;

    class CCA_BT_KEM{
        public:
        CCA_BT_KEM(BTIB* I, Sigma* Sig): I_(I), sigma_(Sig) {}
        CCA_BT_KEM(BTIB* I): I_(I){
            sigma_ = new Sigma();
        }
        CCA_BT_KEM(){
            I_ = new BTIB_KEM1();
            sigma_ = new Sigma();
        }

        BTIBKeyGenOutput KGen(int n, int t, int l, int security_lambda);

        std::pair<Fr, G1> TagKeys();
        proof GenTag(Fr sk_e, const Statement& stmt);
        proof GenTagWithTrapdoor(Fr sk_e, const Statement& stmt);
        bool Verify(const Statement& stmt, const proof& pf);

        std::pair<GT, std::pair<Ciphertext, proof>> Encrypt(const std::shared_ptr<PublicKey>& pk, int j);
        std::pair<GT, Ciphertext> Enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID, Fr r);

        std::pair<bool, std::shared_ptr<idkShare>> PDec(bool b, const SecretKeyComponent& sk_i_b, const Ciphertext& ct);

        bool PVerify(const std::shared_ptr<PublicKey>& pk, const Ciphertext& ciphertext, const std::pair<bool, std::shared_ptr<idkShare>>& d_i);

        GT Comb(const std::vector<int> S,const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, const Ciphertext& ct);

        private:
        static Statement makeStatement(const std::shared_ptr<PublicKey>& pk, const Ciphertext& ct){
            Statement stmt;
            stmt.pk = pk;  // no const_cast, just assign shared_ptr directly
            stmt.ct = ct;
            return stmt;
        }

        BTIB* I_;
        Sigma* sigma_;
    };
}