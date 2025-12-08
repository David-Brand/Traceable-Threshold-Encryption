#include "SignatureScheme.hpp"
#include <array>

namespace proofsystem {

    std::pair<Fr, G1> Sigma::KGen() {
        Fr sk_e;
        sk_e.setByCSPRNG();
        G1 vk_e;
        G1::mul(vk_e, g1(), sk_e);
        return {sk_e, vk_e};
    }

    proof Sigma::signL(Fr randomness, const Statement& stmt) {
        auto const* pk1 = dynamic_cast<const BTIB_KEM1::PublicKey1*>(stmt.pk.get());
        proof pf;

        Fr t0, e1, s1;
        t0.setByCSPRNG();
        e1.setByCSPRNG();
        s1.setByCSPRNG();

        // A0 = Yj^t0, A1 = Zj^t0
        G1::mul(pf.A0, pk1->Y[stmt.ct.j], t0);
        G1::mul(pf.A1, pk1->Z[stmt.ct.j], t0);

        // B = g1^s1 - Q^e1
        G1 g1s1, Qe1;
        G1::mul(g1s1, g1(), s1);
        G1::mul(Qe1, pk1->Q, e1);
        G1::sub(pf.B, g1s1, Qe1);

        Fr e = challenge(stmt, pf.A0, pf.A1, pf.B);
        pf.e1 = e1;
        pf.e0 = e;
        pf.e0 -= pf.e1; // e0 = e - e1

        // s0 = t0 + e0 * randomness  (randomness = r)
        pf.s0 = pf.e0;
        pf.s0 *= randomness;
        pf.s0 += t0;

        pf.s1 = s1;
        return pf;
    }

    proof Sigma::signR(Fr trapdoor, const Statement& stmt) {
        auto const* pk1 = dynamic_cast<const BTIB_KEM1::PublicKey1*>(stmt.pk.get());
        auto const* cc0 = dynamic_cast<const BTIB_KEM1::CiphertextComponent1*>(stmt.ct.c0.get());
        auto const* cc1 = dynamic_cast<const BTIB_KEM1::CiphertextComponent1*>(stmt.ct.c1.get());
        proof pf;

        Fr e0, s0, t1;
        e0.setByCSPRNG();
        s0.setByCSPRNG();
        t1.setByCSPRNG();

        // A0 = Yj^s0 - c0 * e0
        // A1 = Zj^s0 - c1 * e0
        G1 Ys0, Zs0, e0c0, e0c1;
        G1::mul(Ys0,  pk1->Y[stmt.ct.j], s0);
        G1::mul(Zs0,  pk1->Z[stmt.ct.j], s0);
        G1::mul(e0c0, cc0->component, e0);
        G1::mul(e0c1, cc1->component, e0);
        G1::sub(pf.A0, Ys0, e0c0);
        G1::sub(pf.A1, Zs0, e0c1);

        // B = g1^t1
        G1::mul(pf.B, g1(), t1);

        Fr e = challenge(stmt, pf.A0, pf.A1, pf.B);
        pf.e0 = e0;
        pf.e1 = e;
        pf.e1 -= pf.e0; // e1 = e - e0

        pf.s0 = s0; // simulated left response

        // s1 = t1 + e1 * trapdoor  (trapdoor = q)
        pf.s1 = pf.e1;
        pf.s1 *= trapdoor;
        pf.s1 += t1;
        return pf;
    }

    bool Sigma::verify(const Statement& stmt, const proof& pf) {
        auto const* pk1 = dynamic_cast<const BTIB_KEM1::PublicKey1*>(stmt.pk.get());
        auto const* cc0 = dynamic_cast<const BTIB_KEM1::CiphertextComponent1*>(stmt.ct.c0.get());
        auto const* cc1 = dynamic_cast<const BTIB_KEM1::CiphertextComponent1*>(stmt.ct.c1.get());

        Fr e = challenge(stmt, pf.A0, pf.A1, pf.B);

        Fr sum = pf.e0;
        sum += pf.e1;
        if (sum != e) return false;

        G1 lhs, rhs, tmp;

        // left checks:
        G1::mul(lhs, pk1->Y[stmt.ct.j], pf.s0);
        G1::mul(tmp, cc0->component, pf.e0);
        G1::add(rhs, pf.A0, tmp);
        if (lhs != rhs) return false;

        G1::mul(lhs, pk1->Z[stmt.ct.j], pf.s0);
        G1::mul(tmp, cc1->component, pf.e0);
        G1::add(rhs, pf.A1, tmp);
        if (lhs != rhs) return false;

        // right checks:
        G1::mul(lhs, g1(), pf.s1);
        G1::mul(tmp, pk1->Q, pf.e1);
        G1::add(rhs, pf.B, tmp);
        if (lhs != rhs) return false;

        return true;
    }

    Fr Sigma::challenge(const Statement& stmt, const G1& A0, const G1& A1, const G1& B) {
        auto const* pk1 = dynamic_cast<const BTIB_KEM1::PublicKey1*>(stmt.pk.get());
        auto const* cc0 = dynamic_cast<const BTIB_KEM1::CiphertextComponent1*>(stmt.ct.c0.get());
        auto const* cc1 = dynamic_cast<const BTIB_KEM1::CiphertextComponent1*>(stmt.ct.c1.get());
        std::vector<uint8_t> data;
        std::uint64_t v = static_cast<std::uint64_t>(stmt.ct.j);
        for (int i = 0; i < 8; ++i) {
            data.push_back(static_cast<std::uint8_t>(v & 0xff));
            v >>= 8;
        }
        appendVec(data, serializeG1(stmt.ct.ID));
        appendVec(data, serializeG1(pk1->Y[stmt.ct.j]));
        appendVec(data, serializeG1(pk1->Z[stmt.ct.j]));
        appendVec(data, serializeG1(cc0->component));
        appendVec(data, serializeG1(cc1->component));
        appendVec(data, serializeG1(g1()));
        appendVec(data, serializeG1(pk1->Q));
        appendVec(data, serializeG1(A0));
        appendVec(data, serializeG1(A1));
        appendVec(data, serializeG1(B));
        cybozu::Sha256 sha;

        sha.update(data.data(), data.size());

        std::array<uint8_t, 32> out{};
        sha.digest(out.data(), out.size());
        Fr e;
        e.setLittleEndianMod(out.data(), out.size());
        return e;
    }

    std::vector<uint8_t> Sigma::serializeG1(const G1& P) {
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

    std::vector<uint8_t> Sigma::serializeFr(const Fr& x) {
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