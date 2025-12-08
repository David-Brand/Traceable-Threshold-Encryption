#include "CCA_BT_KEM.hpp"

namespace ccakem{

    BTIBKeyGenOutput CCA_BT_KEM::KGen(int n, int t, int l, int security_lambda){
        return I_->keygen(n, t, l, security_lambda);
    }

    std::pair<Fr, G1> CCA_BT_KEM::TagKeys(){
        return sigma_->KGen();
    }

    proof CCA_BT_KEM::GenTag(Fr sk_e, const Statement& stmt){
        return sigma_->signL(sk_e, stmt);
    }
    proof CCA_BT_KEM::GenTagWithTrapdoor(Fr sk_e, const Statement& stmt){
        return sigma_->signR(sk_e, stmt);
    }

    bool CCA_BT_KEM::Verify(const Statement& stmt, const proof& pf){
        return sigma_->verify(stmt, pf);
    }

    std::pair<GT, std::pair<Ciphertext, proof>> CCA_BT_KEM::Encrypt(const std::shared_ptr<PublicKey>& pk, int j){
        auto [sk_e, vk_e] = TagKeys();
        auto [k, ct] = Enc(pk, j, vk_e, sk_e);
        Statement tmp = makeStatement(pk, ct);
        proof pf = GenTag(sk_e, tmp);
        return {k, {ct, pf}};
    }
    std::pair<GT, Ciphertext> CCA_BT_KEM::Enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID, Fr r){
        return I_->enc(pk, j, ID, r);
    }

    std::pair<bool, std::shared_ptr<idkShare>> CCA_BT_KEM::PDec(bool b, const SecretKeyComponent& sk_i_b, const Ciphertext& ct){
        return I_->Didk_j(b, sk_i_b, ct.ID, ct.j);
    }

    bool CCA_BT_KEM::PVerify(const std::shared_ptr<PublicKey>& pk, const Ciphertext& ciphertext, const std::pair<bool, std::shared_ptr<idkShare>>& d_i){
        return I_->PVerify(pk, ciphertext, d_i);
    }

    GT CCA_BT_KEM::Comb(const std::vector<int> S,const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, const Ciphertext& ct){
        auto idk_comb = I_->combIdk_j(S, shares, ct.j);
        return I_->dec(idk_comb, ct);
    }
}