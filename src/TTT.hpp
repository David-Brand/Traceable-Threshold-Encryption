#pragma once

#include "btbf.hpp"
#include "btddh.hpp"
#include "fingerprinting_codes.hpp"
#include <functional>
#include <memory>


namespace ttt{
    using namespace btbf;
    using namespace btddh;
    using namespace fingerprinting;

    struct SecretKey {
        std::vector<SecretKeyComponent> sk_components;
    };

    struct SecretKey_ddh{
        std::vector<SecretKeyComponent_ddh> sk_components;
    };

    struct TTTKeyGenOutput {
        PublicKey* public_key;
        std::vector<SecretKey> parties; // size user_count
        std::vector<std::vector<double>> tracing_key; // tracing key
    };

    struct TTTKeyGenOutput_ddh {
        PublicKey_ddh* public_key;
        std::vector<SecretKey_ddh> parties; // size user_count
        std::vector<std::vector<double>> tracing_key; // tracing key
    };

    class TTT_BTBF {
        public:

        TTT_BTBF(int decryption_threshold = 2, int security_lambda = 4);
        TTTKeyGenOutput keygen(int user_count, int decryption_threshold, int security_lambda, double decoder_error);
        std::pair<SymKey, Ciphertext> enc(const PublicKey* public_key);
        GT dec(const SecretKey& sk, int i, const Ciphertext& c);
        SymKey combine(const std::vector<int>& J, const std::vector<GT>& shares);
        std::vector<std::size_t> trace(const PublicKey* public_key, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext&, const SymKey&)>& D);

        fingerprinting::TardosCodes fingerPrintingCode() { return *fingerPrintingCode_;}

        private:
        fingerprinting::TardosCodes* fingerPrintingCode_;
        int decryption_threshold_, security_lambda_;
        BTBF* BTE;

        int TrD(const PublicKey* public_key, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext&, const SymKey&)>& D);

    };

    class TTT_BTDDH {
        public:

        TTT_BTDDH(int decryption_threshold = 2, int security_lambda = 4);
        TTTKeyGenOutput_ddh keygen(int n, int decryption_threshold, int security_lambda, double decoder_error);
        std::pair<G, Ciphertext_ddh> enc(const PublicKey_ddh* public_key);
        G dec(const SecretKey_ddh& sk, int i, const Ciphertext_ddh& c);
        G combine(const std::vector<int>& J, const std::vector<G>& shares);
        std::vector<std::size_t> trace(const PublicKey_ddh* public_key, std::vector<std::vector<double>> tk, const std::function<bool(const Ciphertext_ddh&, const G&)>& D);
        fingerprinting::TardosCodes fingerPrintingCode() { return *fingerPrintingCode_;}

        private:
        fingerprinting::TardosCodes* fingerPrintingCode_;
        int decryption_threshold_, security_lambda_;
        BTDDH* BTE;

        int TrD(const PublicKey_ddh* public_key, std::size_t j, double N, bool bk, bool b0, bool b1, const std::function<bool(const Ciphertext_ddh&, const G&)>& D);
    };
}
