#include "BTIB.hpp"
#include "Shamir.hpp"

namespace btib {
    BTIBKeyGenOutput BTIB_KEM1::keygen(int n, int t, int l, int security_lambda){
        BTIBKeyGenOutput output;
        auto pk1 = std::make_shared<PublicKey1>();
        pk1->X.resize(l);
        pk1->Y.resize(l);
        pk1->Z.resize(l);

        std::vector<PartySecret> parties(n);
        for(int i=0; i<n; ++i){
            parties[i].left.resize(l);
            parties[i].right.resize(l);
        }

        for(int j=0; j<l; j++){
            Fr alpha, y, z;
            alpha.setByCSPRNG();
            y.setByCSPRNG();
            z.setByCSPRNG();

            // X_j = g1^{alpha*y*z}, Y_j = g1^y, Z_j = g1^z
            G1 X_j, Y_j, Z_j;
            Fr yz; Fr::mul(yz, y, z);
            Fr alpha_yz; Fr::mul(alpha_yz, alpha, yz);
            G1::mul(X_j, g1_, alpha_yz);
            G1::mul(Y_j, g1_, y);
            G1::mul(Z_j, g1_, z);
            pk1->X[j] = X_j;
            pk1->Y[j] = Y_j;
            pk1->Z[j] = Z_j;

            // Shamir shares of alpha at positions 1..n
            auto s = share<Fr>(alpha, n, t);
            for(int i=0; i<n; i++){
                auto left  = std::make_shared<SecretKeyComponent1>();
                auto right = std::make_shared<SecretKeyComponent1>();
                Fr::mul(left->sk,  s[i], z);
                Fr::mul(right->sk, s[i], y);
                parties[i].left[j]  = std::move(left);
                parties[i].right[j] = std::move(right);
            }
        }
        output.pk = pk1;
        output.parties = std::move(parties);
        return output;
    }

    std::pair<GT, Ciphertext> BTIB_KEM1::enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID){
        Fr r;
        r.setByCSPRNG();
        return enc(pk, j, ID, r);
    }
    std::pair<GT, Ciphertext> BTIB_KEM1::enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID, Fr r){
        auto const* pk1 = dynamic_cast<const PublicKey1*>(pk.get());

        // k = e(X_j^r, H(j,ID))
        GT k;
        G1 Xr = pk1->X[j];
        G1::mul(Xr, Xr, r);
        pairing(k, Xr, H(j, ID));

        // c0 = Y_j^r, c1 = Z_j^r
        auto c0 = std::make_shared<CiphertextComponent1>();
        auto c1 = std::make_shared<CiphertextComponent1>();
        G1::mul(c0->component, pk1->Y[j], r);
        G1::mul(c1->component, pk1->Z[j], r);

        Ciphertext ciphertext;
        ciphertext.j = j;
        ciphertext.c0 = std::move(c0);
        ciphertext.c1 = std::move(c1);
        ciphertext.ID = ID;
        return std::make_pair(k, std::move(ciphertext));
    }

    std::pair<bool, std::shared_ptr<idkShare>> BTIB_KEM1::Didk_j(bool b, const SecretKeyComponent& sk_i_b, G1 ID, int j){
        auto const* sk = dynamic_cast<const SecretKeyComponent1*>(&sk_i_b);

        auto sh = std::make_shared<idkShare1>();
        G2::mul(sh->share, H(j, ID), sk->sk);
        return std::make_pair(b, std::move(sh));
    }

    bool BTIB_KEM1::PVerify(const std::shared_ptr<PublicKey>& pk, const Ciphertext& ciphertext, const std::pair<bool, std::shared_ptr<idkShare>>& d_i){
        auto const* pk1 = dynamic_cast<const PublicKey1*>(pk.get());
        auto const* sh = dynamic_cast<const idkShare1*>(d_i.second.get());

        GT left, right;
        pairing(left, g1_, sh->share);
        if (d_i.first) {
            pairing(right, pk1->Z[ciphertext.j], H(ciphertext.j, ciphertext.ID));
        } else {
            pairing(right, pk1->Y[ciphertext.j], H(ciphertext.j, ciphertext.ID));
        }
        return left == right;
    }

    std::shared_ptr<idk> BTIB_KEM1::combIdk_j(const std::vector<int> S, const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, int j){
        auto out = std::make_shared<idk1>();
        out->left.clear();  // ensure zero
        out->right.clear();

        for(size_t k = 0; k < S.size(); ++k){
            int i = S[k];
            auto const* sh = dynamic_cast<const idkShare1*>(shares[k].second.get());

            Fr lambda = lagrangeCoeff(S, i);

            G2 t;
            G2::mul(t, sh->share, lambda);
            if (shares[k].first) {
                G2::add(out->right, out->right, t);
            } else {
                G2::add(out->left, out->left, t);
            }
        }
        return out;
    }

    GT BTIB_KEM1::dec(const std::shared_ptr<idk>& idk_comb, const Ciphertext& ciphertext){
        auto const* acc = dynamic_cast<const idk1*>(idk_comb.get());

        auto const* c0 = dynamic_cast<const CiphertextComponent1*>(ciphertext.c0.get());
        auto const* c1 = dynamic_cast<const CiphertextComponent1*>(ciphertext.c1.get());

        GT lPair, rPair, k;
        pairing(lPair, c0->component, acc->left);
        pairing(rPair, c1->component, acc->right);
        GT::mul(k, lPair, rPair);
        return k;
    }
}