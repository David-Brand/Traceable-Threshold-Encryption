#include "ProofSystems.hpp"
#include <array>

namespace proofsystem {

    std::pair<Zp, G> Sigma::KGen() {
        Zp sk_e;
        sk_e.setByCSPRNG();
        G vk_e;
        G::mul(vk_e, g(), sk_e);
        return {sk_e, vk_e};
    }

    proof Sigma::signL(Zp randomness, const Statement& stmt) {
        auto public_key = stmt.public_key;
        proof pf;

        Zp t0, e1, s1;
        t0.setByCSPRNG();
        e1.setByCSPRNG();
        s1.setByCSPRNG();

        // A0 = Yj^t0, A1 = Zj^t0
        G::mul(pf.A0, public_key->Y[stmt.ct->j], t0);
        G::mul(pf.A1, public_key->Z[stmt.ct->j], t0);

        // B = g^s1 - Q^e1
        G g1s1, Qe1;
        G::mul(g1s1, g(), s1);
        G::mul(Qe1, public_key->Q, e1);
        G::sub(pf.B, g1s1, Qe1);

        Zp e = challenge(stmt, pf.A0, pf.A1, pf.B);
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

    proof Sigma::signR(Zp trapdoor, const Statement& stmt) {
        auto public_key = stmt.public_key;
        auto cc0 = stmt.ct->c0;
        auto cc1 = stmt.ct->c1;
        proof pf;

        Zp e0, s0, t1;
        e0.setByCSPRNG();
        s0.setByCSPRNG();
        t1.setByCSPRNG();

        // A0 = Yj^s0 - c0 * e0
        // A1 = Zj^s0 - c1 * e0
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

        pf.s0 = s0; // simulated left response

        // s1 = t1 + e1 * trapdoor  (trapdoor = q)
        pf.s1 = pf.e1;
        pf.s1 *= trapdoor;
        pf.s1 += t1;
        return pf;
    }

    bool Sigma::verify(const Statement& stmt, const proof& pf) {
        auto public_key = stmt.public_key;
        auto cc0 = stmt.ct->c0;
        auto cc1 = stmt.ct->c1;

        Zp e = challenge(stmt, pf.A0, pf.A1, pf.B);

        Zp sum = pf.e0;
        sum += pf.e1;
        if (sum != e) return false;

        G lhs, rhs, tmp;

        // left checks:
        G::mul(lhs, public_key->Y[stmt.ct->j], pf.s0);
        G::mul(tmp, cc0.comp, pf.e0);
        G::add(rhs, pf.A0, tmp);
        if (lhs != rhs) return false;

        G::mul(lhs, public_key->Z[stmt.ct->j], pf.s0);
        G::mul(tmp, cc1.comp, pf.e0);
        G::add(rhs, pf.A1, tmp);
        if (lhs != rhs) return false;

        // right checks:
        G::mul(lhs, g(), pf.s1);
        G::mul(tmp, public_key->Q, pf.e1);
        G::add(rhs, pf.B, tmp);
        if (lhs != rhs) return false;

        return true;
    }

    Zp Sigma::challenge(const Statement& stmt, const G& A0, const G& A1, const G& B) {
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
        appendVec(data, serializeG1(g()));
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

    std::vector<uint8_t> Sigma::serializeG1(const G& P) {
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

    std::vector<uint8_t> Sigma::serializeFr(const Zp& x) {
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