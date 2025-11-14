#include "btbf.hpp"
#include <mcl/bn.hpp>
#include <openssl/sha.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

using namespace btbf;
using namespace mcl::bls12;

void BTBF::init() {
    // Initialize curve (BLS12-381) and hash-to-curve
    initPairing(mcl::BLS12_381);
    mcl::bn::setMapToMode(MCL_MAP_TO_MODE_HASH_TO_CURVE);

    // Set domain separation tags for hash
    const char* dstG1 = "BTBF:G1";
    const char* dstG2 = "BTBF:G2";
    mcl::bn::setDstG1(dstG1, (int)std::strlen(dstG1));
    mcl::bn::setDstG2(dstG2, (int)std::strlen(dstG2));
}


void BTBF::G1gen(G1& g1) {
    const char* tag = "BTBF:G1:BASE";
    hashAndMapToG1(g1, tag, (size_t)std::strlen(tag));
}
void BTBF::G2gen(G2& g2) {
    const char* tag = "BTBF:G2:BASE";
    hashAndMapToG2(g2, tag, (size_t)std::strlen(tag));
}

//H1 on curve
void BTBF::H1(G1& out, int j) {
    std::string s = "BTBF:H1:" + std::to_string(j);
    hashAndMapToG1(out, s.data(), s.size());
}

//H2 using SHA-256
std::string BTBF::kdf(const GT& x) {
    // Serialize GT to bytes, then SHA-256 -> hex
    std::vector<uint8_t> buf(1024);
    size_t n = x.serialize(buf.data(), buf.size());
    buf.resize(n);

    unsigned char h[SHA256_DIGEST_LENGTH];
    SHA256(buf.data(), buf.size(), h);

    std::ostringstream os;
    for (int i = 0; i < SHA256_DIGEST_LENGTH; i++) {
        os << std::hex << std::setw(2) << std::setfill('0') << (int)h[i];
    }
    return os.str();
}

std::vector<Fr> BTBF::lagrangeAtZero(const std::vector<int>& J) {
    std::vector<Fr> lam(J.size());
    for (size_t a = 0; a < J.size(); ++a) {
        Fr num; num.setStr("1");
        Fr den; den.setStr("1");
        Fr iFr; iFr.setStr(std::to_string(J[a]).c_str());

        for (size_t b = 0; b < J.size(); ++b) {
            if (a == b) continue;
            Fr vFr; vFr.setStr(std::to_string(J[b]).c_str());
            Fr diff = vFr; diff -= iFr;

            num *= vFr;
            den *= diff;
        }
        Fr invDen; Fr::inv(invDen, den);
        lam[a] = num * invDen;
    }
    return lam;
}


void BTBF::KeyGen(int n, int t, int ell, PublicKey& pk, std::vector<PartySecretKey>& sks) {
    Fr alpha, y, z;
    alpha.setByCSPRNG();
    y.setByCSPRNG();
    z.setByCSPRNG();

    // Generators
    G1 g1; G2 g2;
    G1gen(g1);
    G2gen(g2);

    // X = g1^(α*y*z), Y = g1^y, Z = g1^z
    Fr ayz = alpha; ayz *= y; ayz *= z;
    pk.X = g1 * ayz;
    pk.Y = g1 * y;
    pk.Z = g1 * z;

    // Create Shamir sharing of alpha
    std::vector<Fr> coeff(t);
    coeff[0] = alpha;
    for (int k = 1; k < t; ++k) coeff[k].setByCSPRNG();

    auto eval = [&](int x)->Fr {
        Fr X; X.setStr(std::to_string(x).c_str());
        Fr acc; acc.setStr("0");
        Fr pow; pow.setStr("1");
        for (int k = 0; k < t; ++k) {
            acc += coeff[k] * pow;
            pow *= X;
        }
        return acc;
    };

    // Allocate secret keys for parties
    sks.clear();
    sks.resize(n);
    for (int i = 0; i < n; ++i) {
        sks[i].pos.resize(ell);
    }

    //H1(j) for all positions
    std::vector<G1> H1j(ell);
    for (int j = 1; j <= ell; ++j) {
        H1(H1j[j-1], j);
    }

    // For party i, compute s_i and left/right shares at position j
    for (int i = 1; i <= n; ++i) {
        Fr si = eval(i);
        Fr zsi = z * si;
        Fr ysi = y * si;

        for (int j = 1; j <= ell; ++j) {
            Share L;
            L.k0 = H1j[j-1] * zsi; // G1
            L.k1 = g2 * zsi;       // G2

            Share R;
            R.k0 = H1j[j-1] * ysi; // G1
            R.k1 = g2 * ysi;       // G2

            sks[i-1].pos[j-1] = std::make_pair(L, R);
        }
    }
}

void BTBF::Enc(const PublicKey& pk, int j, Ciphertext& ct, SharedKey& k) {
    Fr r, t0, t1;
    r.setByCSPRNG();
    t0.setByCSPRNG();
    t1.setByCSPRNG();

    // Generators and H1(j)
    G2 g2; G2gen(g2);
    G1 Hj; H1(Hj, j);

    GT eg; pairing(eg, pk.X, g2);
    GT W; GT::pow(W, eg, r);
    k.hex = kdf(W);

    ct.c0.u = g2 * t0;
    ct.c0.v = (pk.Y * r) + (Hj * t0);

    ct.c1.u = g2 * t1;
    ct.c1.v = (pk.Z * r) + (Hj * t1);
}

GT BTBF::DecShare(const PartySecretKey& ski, int j, int b, const Ciphertext& ct) {
    const Side& cb = (b == 0) ? ct.c0 : ct.c1;
    const Share& sk = (b == 0) ? ski.pos[j-1].first : ski.pos[j-1].second;

    GT num, den, di;
    pairing(num, cb.v, sk.k1);
    pairing(den, sk.k0, cb.u);
    GT::div(di, num, den);
    return di;
}

SharedKey BTBF::Combine(const std::vector<int>& J, const std::vector<GT>& shares) {
    std::vector<Fr> lam = lagrangeAtZero(J);
    GT W; W.setStr("1");

    for (size_t idx = 0; idx < J.size(); ++idx) {
        GT powi;
        GT::pow(powi, shares[idx], lam[idx]);
        GT::mul(W, W, powi);
    }
    SharedKey k; k.hex = kdf(W);
    return k;
}
