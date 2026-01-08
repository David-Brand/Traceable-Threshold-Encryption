#pragma once

#include "btddh.hpp"
#include <mcl/bls12_381.hpp>
#include <mcl/bn.h>
#include "cybozu/sha2.hpp"

namespace proofsystem {
    using namespace btddh;

    struct proof {
        G A0;
        G A1;
        G B;
        Zp e0;
        Zp e1;
        Zp s0;
        Zp s1;
    };

    struct Statement {
        PublicKey_ddh* public_key;
        Ciphertext_ddh* ct;
    };

    class Sigma {
    public:
        Sigma() = default;

        static std::pair<Zp, G> KGen();

        static proof signL(Zp randomness, const Statement& stmt);
        static proof signR(Zp trapdoor, const Statement& stmt);

        static bool verify(const Statement& stmt, const proof& pf);

        private:
        static Zp challenge(const Statement& stmt, const G& A0, const G& A1, const G& B);
        static G g(){
            G g;
            const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
            g.setStr(g1Str, 16);
            return g;
        }

        static std::vector<uint8_t> serializeG1(const G& P);
        static std::vector<uint8_t> serializeFr(const Zp& x);
        static void appendVec(std::vector<uint8_t>& data, const std::vector<uint8_t>& buf){
            if (buf.empty()) return;
            data.insert(data.end(), buf.begin(), buf.end());
        }
    };
}