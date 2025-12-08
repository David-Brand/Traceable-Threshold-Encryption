#pragma once

#include <mcl/bls12_381.hpp>
#include "BTIB.hpp"
#include "cybozu/sha2.hpp"

namespace proofsystem {
    using namespace btib;

    struct proof {
        G1 A0;
        G1 A1;
        G1 B;
        Fr e0;
        Fr e1;
        Fr s0;
        Fr s1;
    };

    struct Statement {
        std::shared_ptr<PublicKey> pk;
        Ciphertext ct;
    };

    class Sigma {
    public:
        Sigma() = default;

        static std::pair<Fr, G1> KGen();

        static proof signL(Fr randomness, const Statement& stmt);
        static proof signR(Fr trapdoor, const Statement& stmt);

        static bool verify(const Statement& stmt, const proof& pf);

        private:
        static Fr challenge(const Statement& stmt, const G1& A0, const G1& A1, const G1& B);
        static G1 g1(){
            G1 g1;
            const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
            g1.setStr(g1Str, 16);
            return g1;
        }

        static std::vector<uint8_t> serializeG1(const G1& P);
        static std::vector<uint8_t> serializeFr(const Fr& x);
        static void appendVec(std::vector<uint8_t>& data, const std::vector<uint8_t>& buf){
            if (buf.empty()) return;
            data.insert(data.end(), buf.begin(), buf.end());
        }
    };
}