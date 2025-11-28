#ifndef TRACEABLE_TRESHOLD_ENCRYPTION
#define TRACEABLE_TRESHOLD_ENCRYPTION

#include "btbf.hpp"
#include "fingerprinting_codes.hpp"
#include <functional>


namespace TTT{
    using namespace btbf;
    using namespace fingerprinting;

    struct SecretKeyComponent_b {
        BTBF::SecretKeyComponent sk_b;
        bool b; // left (b=0) or right (b=1)
    };

    struct TTTKeyGenOutput {
        BTBF::PublicKey pk;
        std::vector<std::vector<SecretKeyComponent_b>> parties; // size n
    };

    class TTT {
        public:
        using G1 = mcl::bls12::G1;
        using G2 = mcl::bls12::G2;
        using GT = mcl::bls12::GT;
        using Fr = mcl::bls12::Fr;

        TTT();
        TTTKeyGenOutput keygen(int n, int t, int security_lambda, double e) const;
        std::pair<BTBF::SymKey, BTBF::Ciphertext> enc(const BTBF::PublicKey& pk) const;
        GT dec(const SecretKeyComponent_b& sk_b, const BTBF::Ciphertext& c) const;
        BTBF::SymKey combine(const std::vector<int>& J, const std::vector<GT>& shares) const;
        std::size_t trace(BTBF::PublicKey pk, const std::function<bool(const BTBF::Ciphertext&, const BTBF::SymKey&)>& D) const;

        fingerprinting::LogLengthCodes fingerPrintingCode() const { return fingerPrintingCode_;}

        private:
        fingerprinting::LogLengthCodes& fingerPrintingCode_;

        int TrD(BTBF::PublicKey pk, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const BTBF::Ciphertext&, const BTBF::SymKey&)>& D) const;

    };
}

#endif