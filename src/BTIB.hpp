#pragma once

#include <memory>
#include <mcl/bls12_381.hpp>
#include <mcl/curve_type.h>
#include <mcl/bn.h>

namespace btib{
    using G1 = mcl::bls12::G1;
    using G2 = mcl::bls12::G2;
    using GT = mcl::bls12::GT;
    using Fr = mcl::bls12::Fr;


    struct PublicKey{
        G1 Q;
        virtual ~PublicKey() = default;
    };
    struct SecretKeyComponent{
        virtual ~SecretKeyComponent() = default;
    };

    struct PartySecret {
        std::vector<std::shared_ptr<SecretKeyComponent>> left;   // size l
        std::vector<std::shared_ptr<SecretKeyComponent>> right;  // size l
    };

    struct CiphertextComponent{
        virtual ~CiphertextComponent() = default;
    };
    struct idkShare{
        virtual ~idkShare() = default;
    };
    struct idk{
        virtual ~idk() = default;
    };

    struct Ciphertext {
        int j;
        std::shared_ptr<CiphertextComponent> c0;
        std::shared_ptr<CiphertextComponent> c1;
        G1 ID;
    };

    struct BTIBKeyGenOutput {
        std::shared_ptr<PublicKey> pk;
        std::vector<PartySecret> parties; // size n
    };

    class BTIB{
    public:
        BTIB(){
            // Initialize mcl for BLS12-381
            initPairing(mcl::BLS12_381);
            // Use standard hash-to-curve mode
            mcl::setMapToMode(MCL_MAP_TO_MODE_HASH_TO_CURVE);

            // Define g1 and g2 generators
            const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
            g1_.setStr(g1Str, 16);
            const char *g2Str = "1 0x24aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e 0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801 0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be";
            g2_.setStr(g2Str, 16);
        };

        virtual ~BTIB() = default;

        virtual BTIBKeyGenOutput keygen(int n, int t, int l, int security_lambda) = 0;

        virtual std::pair<GT, Ciphertext> enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID) = 0;
        virtual std::pair<GT, Ciphertext> enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID, Fr r) = 0;

        virtual std::pair<bool, std::shared_ptr<idkShare>> Didk_j(bool b, const SecretKeyComponent& sk_i_b, G1 ID, int j) = 0;

        virtual bool PVerify(const std::shared_ptr<PublicKey>& pk, const Ciphertext& ciphertext, const std::pair<bool, std::shared_ptr<idkShare>>& d_i) = 0;

        virtual std::shared_ptr<idk> combIdk_j(const std::vector<int> S, const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, int j) = 0;

        virtual GT dec(const std::shared_ptr<idk>& idk_comb, const Ciphertext& ciphertext) = 0;

        
        G1 g1_; // generator of G1
        G2 g2_; // generator of G2
    private:
    };


    class BTIB_KEM1 : public BTIB{
        public:
        struct PublicKey1 : public PublicKey{
            std::vector<G1> X;
            std::vector<G1> Y;
            std::vector<G1> Z;
        };
        struct SecretKeyComponent1 : public SecretKeyComponent{
            Fr sk;
        };
        struct CiphertextComponent1 : public CiphertextComponent{
            G1 component;
        };
        struct idkShare1 : public idkShare{
            G2 share;
        };
        struct idk1 : public idk{
            G2 left;
            G2 right;
        };

        BTIB_KEM1() : super() {}
        ~BTIB_KEM1(){}

        BTIBKeyGenOutput keygen(int n, int t, int l, int security_lambda) override;

        std::pair<GT, Ciphertext> enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID) override;
        std::pair<GT, Ciphertext> enc(const std::shared_ptr<PublicKey>& pk, int j, G1 ID, Fr r) override;

        std::pair<bool, std::shared_ptr<idkShare>> Didk_j(bool b, const SecretKeyComponent& sk_i_b, G1 ID, int j) override;

        bool PVerify(const std::shared_ptr<PublicKey>& pk, const Ciphertext& ciphertext, const std::pair<bool, std::shared_ptr<idkShare>>& d_i) override;

        std::shared_ptr<idk> combIdk_j(const std::vector<int> S,const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, int j) override;

        GT dec(const std::shared_ptr<idk>& idk_comb, const Ciphertext& ciphertext) override;

        private:
        typedef BTIB super;
        static G2 H(int j, const G1& ID){
            G2 h;
            std::string s =std::to_string(j) + ID.getStr(10);
            hashAndMapToG2(h, s.data(), s.size());
            return h;
        }
        // Lagrange coefficient λ_i^J for points in J
        Fr lagrangeCoeff(const std::vector<int>& J, int i){
            Fr fi = i;
            Fr lambda = 1;

            for (int v : J) {
                if (v == i) continue;

                Fr fv = v;
                Fr sub;
                Fr::sub(sub, fv, fi);
                Fr div;
                Fr::div(div, fv, sub);
                Fr::mul(lambda, lambda, div);
            }
            return lambda;
        }
    };

} // namespace btib
