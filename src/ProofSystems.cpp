// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: MIT

#include "ProofSystems.hpp"
#include <array>

namespace proofsystem {

    proof_Enc_DDH Pi_Enc_DDH::signL(Zp randomness, const Statement_Enc_DDH& stmt) {
        auto public_key = stmt.public_key;
        proof_Enc_DDH pf;

        Zp t0, e1, s1;
        t0.setByCSPRNG();
        e1.setByCSPRNG();
        s1.setByCSPRNG();

        // A0 = Yj^t0, A1 = Zj^t0
        G::mul(pf.A0, public_key->Y[stmt.ct->j], t0);
        G::mul(pf.A1, public_key->Z[stmt.ct->j], t0);

        // B = g^zq - Q^e1
        G g1s1, Qe1;
        G::mul(g1s1, g(), s1);
        G::mul(Qe1, public_key->Q, e1);
        G::sub(pf.B, g1s1, Qe1);

        Zp e = challenge(stmt, pf.A0, pf.A1, pf.B);
        pf.e1 = e1;
        pf.e0 = e;
        pf.e0 -= pf.e1; // e0 = e - e1

        // zr = t0 + e0 * randomness  (randomness = r)
        pf.zr = pf.e0;
        pf.zr *= randomness;
        pf.zr += t0;

        pf.zq = s1;
        return pf;
    }

    proof_Enc_DDH Pi_Enc_DDH::signR(Zp trapdoor, const Statement_Enc_DDH& stmt) {
        auto public_key = stmt.public_key;
        auto cc0 = stmt.ct->c0;
        auto cc1 = stmt.ct->c1;
        proof_Enc_DDH pf;

        Zp e0, s0, t1;
        e0.setByCSPRNG();
        s0.setByCSPRNG();
        t1.setByCSPRNG();

        // A0 = Yj^zr - c0 * e0
        // A1 = Zj^zr - c1 * e0
        G Ys0, Zs0, e0c0, e0c1;
        G::mul(Ys0,  public_key->Y[stmt.ct->j], s0);
        G::mul(Zs0,  public_key->Z[stmt.ct->j], s0);
        G::mul(e0c0, cc0.comp, e0);
        G::mul(e0c1, cc1.comp, e0);
        G::sub(pf.A0, Ys0, e0c0);
        G::sub(pf.A1, Zs0, e0c1);

        // B = g^t1
        G::mul(pf.B, g(), t1);
        Zp e = challenge(stmt, pf.A0, pf.A1, pf.B);
        pf.e0 = e0;
        pf.e1 = e;
        pf.e1 -= pf.e0; // e1 = e - e0

        pf.zr = s0; // simulated left response

        // zq = t1 + e1 * trapdoor  (trapdoor = q)
        pf.zq = pf.e1;
        pf.zq *= trapdoor;
        pf.zq += t1;
        return pf;
    }

    bool Pi_Enc_DDH::verify(const Statement_Enc_DDH& stmt, const proof_Enc_DDH& pf) {
        auto public_key = stmt.public_key;
        auto cc0 = stmt.ct->c0;
        auto cc1 = stmt.ct->c1;

        Zp e = challenge(stmt, pf.A0, pf.A1, pf.B);

        Zp sum = pf.e0;
        sum += pf.e1;
        if (sum != e) return false;

        G lhs, rhs, tmp;

        // left checks:
        G::mul(lhs, public_key->Y[stmt.ct->j], pf.zr);
        G::mul(tmp, cc0.comp, pf.e0);
        G::add(rhs, pf.A0, tmp);
        if (lhs != rhs) return false;

        G::mul(lhs, public_key->Z[stmt.ct->j], pf.zr);
        G::mul(tmp, cc1.comp, pf.e0);
        G::add(rhs, pf.A1, tmp);
        if (lhs != rhs) return false;

        // right checks:
        G::mul(lhs, g(), pf.zq);
        G::mul(tmp, public_key->Q, pf.e1);
        G::add(rhs, pf.B, tmp);
        if (lhs != rhs) return false;

        return true;
    }

    Zp Pi_Enc_DDH::challenge(const Statement_Enc_DDH& stmt, const G& A0, const G& A1, const G& B) {
        auto public_key = stmt.public_key;
        auto cc0 = stmt.ct->c0;
        auto cc1 = stmt.ct->c1;
        std::vector<uint8_t> data;
        std::uint64_t v = static_cast<std::uint64_t>(stmt.ct->j);
        for (int i = 0; i < 8; ++i) {
            data.push_back(static_cast<std::uint8_t>(v & 0xff));
            v >>= 8;
        }
        appendVec(data, serializeG1(public_key->Y[stmt.ct->j]));
        appendVec(data, serializeG1(public_key->Z[stmt.ct->j]));
        appendVec(data, serializeG1(cc0.comp));
        appendVec(data, serializeG1(cc1.comp));
        appendVec(data, serializeG1(public_key->Q));
        appendVec(data, serializeG1(A0));
        appendVec(data, serializeG1(A1));
        appendVec(data, serializeG1(B));
        cybozu::Sha256 sha;

        sha.update(data.data(), data.size());

        std::array<uint8_t, 32> out{};
        sha.digest(out.data(), out.size());
        Zp e;
        e.setLittleEndianMod(out.data(), out.size());
        return e;
    }

    std::vector<uint8_t> Pi_Enc_DDH::serializeG1(const G& P) {
        const size_t size = mclBn_getG1ByteSize();
        std::vector<uint8_t> out(size);
        if (P.isZero()) {
            std::fill(out.begin(), out.end(), 0);
            return out;
        }
        G Q;
        G::normalize(Q, P);

        const size_t fpSize = mclBn_getFpByteSize();
        size_t written = Q.x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = P.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        std::vector<uint8_t> yBuf(fpSize);
        written = Q.y.serialize(yBuf.data(), yBuf.size());
        if (written == 0) {
            throw std::runtime_error("Fp::serialize for G1.y failed");
        }
        bool yIsOdd = (yBuf[0] & 1) != 0;
        if (yIsOdd) {
            out[fpSize - 1] |= 0x80;
        }
        return out;
    }

    std::vector<uint8_t> Pi_Enc_DDH::serializeFr(const Zp& x) {
        const size_t size = mclBn_getFrByteSize();
        std::vector<uint8_t> out(size);
        size_t written = x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = x.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        return out;
    }

    proof_Dec_DDH Pi_Dec_DDH::signL(Zp witness, const Statement_Dec_DDH& stmt) {
        proof_Dec_DDH pf;

        Zp sl, er, zr;
        sl.setByCSPRNG();
        er.setByCSPRNG();
        zr.setByCSPRNG();

        // A_l=c_0^s_l
        G::mul(pf.Al, stmt.ct->c0.comp, sl);

        // A_r = c_1^z_r / d_i^e_r
        G c1zr, die1;
        G::mul(c1zr, stmt.ct->c1.comp, zr);
        G::mul(die1, stmt.decryption_share, er);
        G::sub(pf.Ar, c1zr, die1);

        Zp e = challenge(stmt, pf.Al, pf.Ar);
        pf.er = er;
        pf.el = e;
        pf.el -= pf.er; // el = e - er

        // zl = sl + el * witness
        pf.zl = pf.el;
        pf.zl *= witness;
        pf.zl += sl;

        pf.zr = zr;
        return pf;
    }

    proof_Dec_DDH Pi_Dec_DDH::signR(Zp witness, const Statement_Dec_DDH& stmt) {
        proof_Dec_DDH pf;

        Zp sr, el, zl;
        sr.setByCSPRNG();
        el.setByCSPRNG();
        zl.setByCSPRNG();

        // A_r=c_1^s_r
        G::mul(pf.Ar, stmt.ct->c1.comp, sr);

        // A_l = c_0^z_l / d_i^e_l
        G c0zl, die1;
        G::mul(c0zl, stmt.ct->c0.comp, zl);
        G::mul(die1, stmt.decryption_share, el);
        G::sub(pf.Al, c0zl, die1);

        Zp e = challenge(stmt, pf.Al, pf.Ar);
        pf.el = el;
        pf.er = e;
        pf.er -= pf.el; // er = e - el

        // zr = sr + er * witness
        pf.zr = pf.er;
        pf.zr *= witness;
        pf.zr += sr;

        pf.zl = zl;
        return pf;
    }

    bool Pi_Dec_DDH::verify(const Statement_Dec_DDH& stmt, const proof_Dec_DDH& pf) {
        Zp e = challenge(stmt, pf.Al, pf.Ar);

        Zp sum = pf.el;
        sum += pf.er;
        if (sum != e) return false;

        G lhs, rhs, tmp;

        // left check:
        G::mul(lhs, stmt.ct->c0.comp, pf.zl);
        G::mul(tmp, stmt.decryption_share, pf.el);
        G::add(rhs, pf.Al, tmp);
        if (lhs != rhs) return false;

        // right check:
        G::mul(lhs, stmt.ct->c1.comp, pf.zr);
        G::mul(tmp, stmt.decryption_share, pf.er);
        G::add(rhs, pf.Ar, tmp);
        if (lhs != rhs) return false;

        return true;
    }

    Zp Pi_Dec_DDH::challenge(const Statement_Dec_DDH& stmt, const G& Al, const G& Ar) {
        std::vector<uint8_t> data;
        std::uint64_t v = static_cast<std::uint64_t>(stmt.ct->j);
        for (int i = 0; i < 8; ++i) {
            data.push_back(static_cast<std::uint8_t>(v & 0xff));
            v >>= 8;
        }
        appendVec(data, serializeG1(stmt.decryption_share));
        appendVec(data, serializeG1(stmt.ct->c0.comp));
        appendVec(data, serializeG1(stmt.ct->c1.comp));
        appendVec(data, serializeG1(Al));
        appendVec(data, serializeG1(Ar));
        cybozu::Sha256 sha;

        sha.update(data.data(), data.size());

        std::array<uint8_t, 32> out{};
        sha.digest(out.data(), out.size());
        Zp e;
        e.setLittleEndianMod(out.data(), out.size());
        return e;
    }

    std::vector<uint8_t> Pi_Dec_DDH::serializeG1(const G& P) {
        const size_t size = mclBn_getG1ByteSize();
        std::vector<uint8_t> out(size);
        if (P.isZero()) {
            std::fill(out.begin(), out.end(), 0);
            return out;
        }
        G Q;
        G::normalize(Q, P);

        const size_t fpSize = mclBn_getFpByteSize();
        size_t written = Q.x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = P.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        std::vector<uint8_t> yBuf(fpSize);
        written = Q.y.serialize(yBuf.data(), yBuf.size());
        if (written == 0) {
            throw std::runtime_error("Fp::serialize for G1.y failed");
        }
        bool yIsOdd = (yBuf[0] & 1) != 0;
        if (yIsOdd) {
            out[fpSize - 1] |= 0x80;
        }
        return out;
    }

    std::vector<uint8_t> Pi_Dec_DDH::serializeFr(const Zp& x) {
        const size_t size = mclBn_getFrByteSize();
        std::vector<uint8_t> out(size);
        size_t written = x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = x.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        return out;
    }

    proof_Enc_BTBF Pi_Enc_BTBF::signL(Zp randomness, Zp t_0, Zp t_1, const Statement_Enc_BTBF& stmt) {
        auto public_key = stmt.public_key;
        proof_Enc_BTBF pf;

        Zp sr, s0, s1, e1, zq;
        sr.setByCSPRNG();
        s0.setByCSPRNG();
        s1.setByCSPRNG();
        e1.setByCSPRNG();
        zq.setByCSPRNG();

        // Au0 = g2^s0, Au1 = g2^s1
        G2::mul(pf.Au0, g2(), s0);
        G2::mul(pf.Au1, g2(), s1);

        // Av0 = Y^sr * H1(j) ^ s0, Av1 = Z^sr * H1(j) ^ s1
        G1 Hj = BTBF::H1(stmt.ct->j);
        G1 tmp;
        G1::mul(tmp, Hj, s0);
        G1::mul(pf.Av0, public_key->Y, sr);
        G1::add(pf.Av0, pf.Av0, tmp);

        G1::mul(tmp, Hj, s1);
        G1::mul(pf.Av1, public_key->Z, sr);
        G1::add(pf.Av1, pf.Av1, tmp);

        // AQ = g1^zq / Q^e1
        G g1zq, Qe1;
        G::mul(g1zq, g1(), zq);
        G::mul(Qe1, public_key->Q, e1);
        G::sub(pf.AQ, g1zq, Qe1);

        Zp e = challenge(stmt, pf.Au0, pf.Au1, pf.Av0, pf.Av1, pf.AQ);
        pf.e1 = e1;
        pf.e0 = e;
        pf.e0 -= pf.e1; // e0 = e - e1

        // zr = sr + e0 * randomness
        pf.zr = pf.e0;
        pf.zr *= randomness;
        pf.zr += sr;

        // z0 = s0 + e0 * t_0
        pf.z0 = pf.e0;
        pf.z0 *= t_0;
        pf.z0 += s0;

        // z1 = s1 + e0 * t_1
        pf.z1 = pf.e0;
        pf.z1 *= t_1;
        pf.z1 += s1;

        pf.zq = zq;
        return pf;
    }

    proof_Enc_BTBF Pi_Enc_BTBF::signR(Zp trapdoor, const Statement_Enc_BTBF& stmt) {
        auto public_key = stmt.public_key;
        proof_Enc_BTBF pf;

        Zp e0, zr, z0, z1, sq;
        e0.setByCSPRNG();
        zr.setByCSPRNG();
        z0.setByCSPRNG();
        z1.setByCSPRNG();
        sq.setByCSPRNG();

        // AQ = g1^sq
        G::mul(pf.AQ, g1(), sq);

        // Au0 = g2^z0 / u0 ^ e0
        // Au1 = g2^z1 / u1 ^ e0
        G2 g2z0, g2z1, e0u0, e0u1;
        G2::mul(g2z0, g2(), z0);
        G2::mul(g2z1, g2(), z1);
        G2::mul(e0u0, stmt.ct->c0.u, e0);
        G2::mul(e0u1, stmt.ct->c1.u, e0);
        G2::sub(pf.Au0, g2z0, e0u0);
        G2::sub(pf.Au1, g2z1, e0u1);

        // Av0 = Y^zr * H1(j) ^ z0 / v0 ^ e0
        // Av1 = Z^zr * H1(j) ^ z1 / v1 ^ e0
        G1 Hj = BTBF::H1(stmt.ct->j);
        G1 Yzr, Zzr, Hjz0, Hjz1, e0v0, e0v1;
        G1::mul(Yzr, public_key->Y, zr);
        G1::mul(Zzr, public_key->Z, zr);
        G1::mul(Hjz0, Hj, z0);
        G1::mul(Hjz1, Hj, z1);
        G1::mul(e0v0, stmt.ct->c0.v, e0);
        G1::mul(e0v1, stmt.ct->c1.v, e0);
        G1::add(pf.Av0, Yzr, Hjz0);
        G1::sub(pf.Av0, pf.Av0, e0v0);
        G1::add(pf.Av1, Zzr, Hjz1);
        G1::sub(pf.Av1, pf.Av1, e0v1);

        Zp e = challenge(stmt, pf.Au0, pf.Au1, pf.Av0, pf.Av1, pf.AQ);
        pf.e0 = e0;
        pf.e1 = e;
        pf.e1 -= pf.e0; // e1 = e - e0

        pf.zr = zr; // simulated left response
        pf.z0 = z0;
        pf.z1 = z1;

        // zq = sq + e1 * trapdoor  (trapdoor = q)
        pf.zq = pf.e1;
        pf.zq *= trapdoor;
        pf.zq += sq;
        return pf;
    }

    bool Pi_Enc_BTBF::verify(const Statement_Enc_BTBF& stmt, const proof_Enc_BTBF& pf) {
        auto public_key = stmt.public_key;

        Zp e = challenge(stmt, pf.Au0, pf.Au1, pf.Av0, pf.Av1, pf.AQ);

        Zp sum = pf.e0;
        sum += pf.e1;
        if (sum != e) return false;

        G2 lhs2, rhs2, tmp2;
        G1 lhs, rhs, tmpl1, tmpl2, tmpr;
        G1 Hj = BTBF::H1(stmt.ct->j);

        // left checks:
        G2::mul(lhs2, g2(), pf.z0);
        G2::mul(tmp2, stmt.ct->c0.u, pf.e0);
        G2::add(rhs2, pf.Au0, tmp2);
        if (lhs2 != rhs2) return false;

        G2::mul(lhs2, g2(), pf.z1);
        G2::mul(tmp2, stmt.ct->c1.u, pf.e0);
        G2::add(rhs2, pf.Au1, tmp2);
        if (lhs2 != rhs2) return false;

        G1::mul(tmpl1, public_key->Y, pf.zr);
        G1::mul(tmpl2, Hj, pf.z0);
        G1::add(lhs, tmpl1, tmpl2);
        G1::mul(tmpr, stmt.ct->c0.v, pf.e0);
        G1::add(rhs, pf.Av0, tmpr);
        if (lhs != rhs) return false;

        G1::mul(tmpl1, public_key->Z, pf.zr);
        G1::mul(tmpl2, Hj, pf.z1);
        G1::add(lhs, tmpl1, tmpl2);
        G1::mul(tmpr, stmt.ct->c1.v, pf.e0);
        G1::add(rhs, pf.Av1, tmpr);
        if (lhs != rhs) return false;

        // right checks:
        G1::mul(lhs, g1(), pf.zq);
        G1::mul(tmpr, public_key->Q, pf.e1);
        G1::add(rhs, pf.AQ, tmpr);
        if (lhs != rhs) return false;

        return true;
    }

    Zp Pi_Enc_BTBF::challenge(const Statement_Enc_BTBF& stmt, const G2& Au0, const G2& Au1, const G1& Av0, const G1& Av1, const G1& AQ) {
        auto public_key = stmt.public_key;

        std::vector<uint8_t> data;
        std::uint64_t v = static_cast<std::uint64_t>(stmt.ct->j);
        for (int i = 0; i < 8; ++i) {
            data.push_back(static_cast<std::uint8_t>(v & 0xff));
            v >>= 8;
        }
        appendVec(data, serializeG1(public_key->Y));
        appendVec(data, serializeG1(public_key->Z));
        appendVec(data, serializeG1(public_key->X));
        appendVec(data, serializeG1(public_key->Q));
        appendVec(data, serializeG2(stmt.ct->c0.u));
        appendVec(data, serializeG2(stmt.ct->c1.u));
        appendVec(data, serializeG1(stmt.ct->c0.v));
        appendVec(data, serializeG1(stmt.ct->c1.v));
        appendVec(data, serializeG2(Au0));
        appendVec(data, serializeG2(Au1));
        appendVec(data, serializeG1(Av0));
        appendVec(data, serializeG1(Av1));
        appendVec(data, serializeG1(AQ));
        cybozu::Sha256 sha;

        sha.update(data.data(), data.size());

        std::array<uint8_t, 32> out{};
        sha.digest(out.data(), out.size());
        Zp e;
        e.setLittleEndianMod(out.data(), out.size());
        return e;
    }

    std::vector<uint8_t> Pi_Enc_BTBF::serializeG1(const G1& P) {
        const size_t size = mclBn_getG1ByteSize();
        std::vector<uint8_t> out(size);
        if (P.isZero()) {
            std::fill(out.begin(), out.end(), 0);
            return out;
        }
        G1 Q;
        G1::normalize(Q, P);

        const size_t fpSize = mclBn_getFpByteSize();
        size_t written = Q.x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = P.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        std::vector<uint8_t> yBuf(fpSize);
        written = Q.y.serialize(yBuf.data(), yBuf.size());
        if (written == 0) {
            throw std::runtime_error("Fp::serialize for G1.y failed");
        }
        bool yIsOdd = (yBuf[0] & 1) != 0;
        if (yIsOdd) {
            out[fpSize - 1] |= 0x80;
        }
        return out;
    }

    std::vector<uint8_t> Pi_Enc_BTBF::serializeG2(const G2& P) {
        const size_t size = mclBn_getG2ByteSize();
        std::vector<uint8_t> out(size);
        if (P.isZero()) {
            std::fill(out.begin(), out.end(), 0);
            return out;
        }
        G2 Q;
        G2::normalize(Q, P);

        const size_t fpSize = mclBn_getFpByteSize();
        size_t written = Q.x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = P.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        std::vector<uint8_t> yBuf(fpSize);
        written = Q.y.serialize(yBuf.data(), yBuf.size());
        if (written == 0) {
            throw std::runtime_error("Fp::serialize for G1.y failed");
        }
        bool yIsOdd = (yBuf[0] & 1) != 0;
        if (yIsOdd) {
            out[fpSize - 1] |= 0x80;
        }
        return out;
    }

    std::vector<uint8_t> Pi_Enc_BTBF::serializeFr(const Zp& x) {
        const size_t size = mclBn_getFrByteSize();
        std::vector<uint8_t> out(size);
        size_t written = x.serialize(out.data(), out.size());
        if (written == 0) {
            std::string s = x.getStr(10);
            return std::vector<uint8_t>(s.begin(), s.end());
        }
        return out;
    }

} // namespace proofsystem