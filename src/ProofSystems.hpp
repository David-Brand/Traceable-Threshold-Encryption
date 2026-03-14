#pragma once

// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: GPL-3.0-or-later

#include "btddh.hpp"
#include "btbf.hpp"
#include <mcl/bls12_381.hpp>
#include <mcl/bn.h>
#include "cybozu/sha2.hpp"

using namespace btddh;
using namespace btbf;

namespace proofsystem {

    struct proof_Enc_DDH {
        G A0;
        G A1;
        G B;
        Zp e0;
        Zp e1;
        Zp zr;
        Zp zq;
    };

    struct Statement_Enc_DDH {
        PublicKey_ddh* public_key;
        Ciphertext_ddh* ct;
    };

    struct proof_Dec_DDH {
        G Al;
        G Ar;
        Zp el;
        Zp er;
        Zp zl;
        Zp zr;
    };

    struct Statement_Dec_DDH {
        G decryption_share;
        Ciphertext_ddh* ct;
    };

    struct proof_Enc_BTBF {
        G2 Au0;
        G2 Au1;
        G1 Av0;
        G1 Av1;
        G1 AQ;
        Zp e0;
        Zp e1;
        Zp zr;
        Zp z0;
        Zp z1;
        Zp zq;
    };

    struct Statement_Enc_BTBF {
        PublicKey* public_key;
        Ciphertext* ct;
    };

    class Pi_Enc_DDH {
    public:
        Pi_Enc_DDH() = default;

        static proof_Enc_DDH signL(Zp randomness, const Statement_Enc_DDH& stmt);
        static proof_Enc_DDH signR(Zp trapdoor, const Statement_Enc_DDH& stmt);

        static bool verify(const Statement_Enc_DDH& stmt, const proof_Enc_DDH& pf);

        private:
        static Zp challenge(const Statement_Enc_DDH& stmt, const G& A0, const G& A1, const G& B);
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

    class Pi_Dec_DDH {
    public:
        Pi_Dec_DDH() = default;

        static proof_Dec_DDH signL(Zp randomness, const Statement_Dec_DDH& stmt);
        static proof_Dec_DDH signR(Zp trapdoor, const Statement_Dec_DDH& stmt);

        static bool verify(const Statement_Dec_DDH& stmt, const proof_Dec_DDH& pf);

        private:
        static Zp challenge(const Statement_Dec_DDH& stmt, const G& Al, const G& Ar);
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

    class Pi_Enc_BTBF {
    public:
        Pi_Enc_BTBF() = default;

        static proof_Enc_BTBF signL(Zp randomness, Zp t_0, Zp t_1, const Statement_Enc_BTBF& stmt);
        static proof_Enc_BTBF signR(Zp trapdoor, const Statement_Enc_BTBF& stmt);

        static bool verify(const Statement_Enc_BTBF& stmt, const proof_Enc_BTBF& pf);

        private:
        static Zp challenge(const Statement_Enc_BTBF& stmt, const G2& Au0, const G2& Au1, const G1& Av0, const G1& Av1, const G1& AQ);
        static G1 g1(){
            G1 g;
            const char *g1Str = "1 0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb 0x08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1";
            g.setStr(g1Str, 16);
            return g;
        }
        static G2 g2(){
            G2 g;
            const char *g2Str = "1 0x24aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e 0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801 0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be";
            g.setStr(g2Str, 16);
            return g;
        }

        static std::vector<uint8_t> serializeG1(const G1& P);
        static std::vector<uint8_t> serializeG2(const G2& P);
        static std::vector<uint8_t> serializeFr(const Zp& x);
        static void appendVec(std::vector<uint8_t>& data, const std::vector<uint8_t>& buf){
            if (buf.empty()) return;
            data.insert(data.end(), buf.begin(), buf.end());
        }
    };
}
