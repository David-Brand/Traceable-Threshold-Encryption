#pragma once

#include "CCA_BT_KEM.hpp"
#include "fingerprinting_codes.hpp"
#include <functional>


namespace ttt{
    using namespace ccakem;
    using namespace fingerprinting;

    struct SecretKey {
        std::vector<std::shared_ptr<SecretKeyComponent>> sk_components;
    };

    struct TTTKeyGenOutput {
        std::shared_ptr<PublicKey> pk;
        std::vector<SecretKey> parties; // size n
        std::vector<std::vector<double>> tk; // tracing key
    };

    class TTT {
        public:

        TTT(int t = 2, int security_lambda = 4);
        TTTKeyGenOutput keygen(int n, int t, int security_lambda, double e);
        std::pair<GT, std::pair<Ciphertext, proof>> enc(const std::shared_ptr<PublicKey>& pk);
        std::pair<bool, std::shared_ptr<idkShare>> dec(const SecretKey& sk, int i, const Ciphertext& c);
        GT combine(const std::vector<int>& J, const std::vector<std::pair<bool, std::shared_ptr<idkShare>>>& shares, const Ciphertext& c);
        std::vector<std::size_t> trace(const std::shared_ptr<PublicKey>& pk, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext&, const GT&)>& D);

        fingerprinting::TardosCodes fingerPrintingCode() { return *fingerPrintingCode_;}

        private:
        fingerprinting::TardosCodes* fingerPrintingCode_;
        int t_, security_lambda_;
        CCA_BT_KEM* BTE;

        int TrD(const std::shared_ptr<PublicKey>& pk, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext&, const GT&)>& D);

    };
}
