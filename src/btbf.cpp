#include "btbf.hpp"

#include <algorithm>
#include <stdexcept>
#include <iostream> // optional, for debugging/logging

using namespace mcl::bn;

namespace btbf {

bool BTBF::initialized_ = false;
BTBF::G1 BTBF::g1_;
BTBF::G2 BTBF::g2_;

void BTBF::init(){
    if (initialized_) return;
    // Initialize mcl for BLS12-381
    initPairing(mcl::BLS12_381);
    // Use standard hash-to-curve mode
    setMapToMode(MCL_MAP_TO_MODE_HASH_TO_CURVE);

    // Define g1 and g2 generators
    g1_ = getG1basePoint();
    const char *g2Str = "1 0x24aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 0x13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e 0x0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801 0x0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be";
    g2_.setStr(g2Str, 16);

    initialized_ = true;
}

// H1(j) : N -> G1
BTBF::G1 BTBF::H1(int j){
    G1 P;
    std::string msg = std::to_string(j);
    hashAndMapToG1(P, msg.data(), msg.size());
    return P;
}

// H2(W) : GT -> {0,1}^256
BTBF::SymKey BTBF::H2(const GT& w){
    // Serialize GT element
    uint8_t buf[512];
    size_t n = w.serialize(buf, sizeof(buf));

    // Hash to Fr
    Fr x;
    x.setHashOf(buf, n);

    SymKey key{};
    x.getLittleEndian(key.data(), key.size());
    return key;
}

// Lagrange coefficient λ_i^J at 0 for points in J
BTBF::Fr BTBF::lagrangeCoeffAtZero(const std::vector<int>& J, int i){
    Fr num = 1;
    Fr den = 1;

    Fr fi = i;

    for (int v : J) {
        if (v == i) continue;

        Fr fv = v;
        // num = num * fv;
        Fr::mul(num, num, fv);
        // tmp = fv - fi;
        Fr tmp;
        Fr::sub(tmp, fv, fi);
        // den = den * tmp;
        Fr::mul(den, den, tmp);
    }
    Fr denInv;
    Fr::inv(denInv, den);

    Fr lambda;
    // lambda = num * denInv;
    Fr::mul(lambda, num, denInv);

    return lambda;
}

// BTBF.KeyGen(1^λ, n, t, ℓ)
BTBF::KeyGenOutput BTBF::keygen(int n, int t, int ell)
{
    if (!initialized_) {
        throw std::runtime_error("BTBF::init() was not called");
    }
    if (t <= 0 || t > n) {
        throw std::invalid_argument("t must be in [1, n]");
    }
    if (ell <= 0) {
        throw std::invalid_argument("ell must be > 0");
    }

    KeyGenOutput out;
    out.n = n;
    out.t = t;
    out.ell = ell;
    out.parties.resize(n);

    // 1. Sample α, y, z in Z_p (we use Fr as scalar field)
    Fr alpha, y, z;
    alpha.setByCSPRNG();
    y.setByCSPRNG();
    z.setByCSPRNG();

    // 1. Set X ← g1^{α y z}, Y ← g1^y, Z ← g1^z.
    Fr yz;
    Fr::mul(yz, y, z);

    Fr alpha_yz;
    Fr::mul(alpha_yz, alpha, yz);

    G1 X;
    G1::mul(X, g1_, alpha_yz);

    G1 Y;
    G1::mul(Y, g1_, y);

    G1 Z;
    G1::mul(Z, g1_, z);

    out.pk.X = X;
    out.pk.Y = Y;
    out.pk.Z = Z;

    // Shamir t-out-of-n secret sharing s_1,...,s_n of α
    // polynomial f(x) = alpha + a1 x + ... + a_{t-1} x^{t-1}
    std::vector<Fr> coef(t);
    coef[0] = alpha;
    for (int k = 1; k < t; k++) {
        coef[k].setByCSPRNG();
    }

    std::vector<Fr> s(n);
    for (int i = 0; i < n; i++) {
        int idx = i + 1;    // shares at positions 1..n
        Fr x = idx;

        Fr val = coef[0];
        Fr pow = x;         // x^1

        for (int k = 1; k < t; k++) {
            Fr term;
            // term = coef[k] * pow;
            Fr::mul(term, coef[k], pow);

            // val = val + term;
            Fr::add(val, val, term);

            // pow = pow * x;
            Fr::mul(pow, pow, x);
        }
        s[i] = val; // s_{i+1}
    }

    // Precompute H1(j) for j=1..ell
    std::vector<G1> h1(ell);
    for (int j = 0; j < ell; j++) {
        h1[j] = H1(j + 1);
    }

    // 4. For i=1..n and j=1..ℓ:
    //    k(0)_0 ← H1(j)^{z s_i}, k(0)_1 ← g2^{z s_i}, sk_i,0^{(j)} = (k(0)_0, k(0)_1)
    //    k(1)_0 ← H1(j)^{y s_i}, k(1)_1 ← g2^{y s_i}, sk_i,1^{(j)} = (k(1)_0, k(1)_1)
    for (int i = 0; i < n; i++) {
        PartySecret ps;
        ps.left.resize(ell);
        ps.right.resize(ell);

        Fr zsi;
        Fr::mul(zsi, z, s[i]);

        Fr ysi;
        Fr::mul(ysi, y, s[i]);

        for (int j = 0; j < ell; j++) {
            // left key
            SecretKeyComponent left;
            left.k0 = h1[j];
            // left.k0 = H1(j)^{z s_i}
            G1::mul(left.k0, left.k0, zsi);
            // left.k1 = g2^{z s_i}
            G2::mul(left.k1, g2_, zsi);
            ps.left[j] = left;

            // right key
            SecretKeyComponent right;
            right.k0 = h1[j];
            // right.k0 = H1(j)^{y s_i}
            G1::mul(right.k0, right.k0, ysi);
            // right.k1 = g2^{y s_i}
            G2::mul(right.k1, g2_, ysi);
            ps.right[j] = right;
        }

        out.parties[i] = std::move(ps);
    }
    //pkc is omitted (as it is ⊥)
    return out;
}

// BTBF.Enc(pk, j)
std::pair<BTBF::SymKey, BTBF::Ciphertext> BTBF::encaps(const PublicKey& pk, int j){
    if (!initialized_) {
        throw std::runtime_error("BTBF::init() was not called");
    }
    if (j <= 0) {
        throw std::invalid_argument("j must be >= 1");
    }
    // Sample r, t0, t1
    Fr r, t0, t1;
    r.setByCSPRNG();
    t0.setByCSPRNG();
    t1.setByCSPRNG();

    // W = e(X, g2)^r
    GT W;
    pairing(W, pk.X, g2_);
    GT::pow(W, W, r);
    // k = H2(W)
    SymKey k = H2(W);

    // H1(j)
    G1 Hj = H1(j);

    // c0 : u0 = g2^{t0}, v0 = Y^r * Hj^{t0}
    CiphertextComponent c0;
    // c0.u = g2_^t0
    G2::mul(c0.u, g2_, t0);

    // Y_r = Y^r
    G1 Y_r;
    G1::mul(Y_r, pk.Y, r);
    // Hj_t0 = Hj^{t0}
    G1 Hj_t0;
    G1::mul(Hj_t0, Hj, t0);
    // c0.v = Y_r * Hj_t0
    G1::add(c0.v, Y_r, Hj_t0);

    // c1 : u1 = g2^{t1}, v1 = Z^r * Hj^{t1}
    CiphertextComponent c1;
    // c1.u = g2_^t1
    G2::mul(c1.u, g2_, t1);

    // Z_r = Z^r
    G1 Z_r;
    G1::mul(Z_r, pk.Z, r);
    // Hj_t1 = Hj^{t1}
    G1 Hj_t1;
    G1::mul(Hj_t1, Hj, t1);
    // c1.v = Z_r * Hj_t1
    G1::add(c1.v, Z_r, Hj_t1);

    Ciphertext c;
    c.c0 = c0;
    c.c1 = c1;
    c.j  = j;
    return { k, c };
}

// BTBF.Dec(j, sk_i,b^{(j)}, c_b)
BTBF::GT BTBF::decShare(const SecretKeyComponent& sk_i_b, const CiphertextComponent& c_b){
    GT e1, e2;
    pairing(e1, c_b.v, sk_i_b.k1);  // e(v, k1)
    pairing(e2, sk_i_b.k0, c_b.u);  // e(k0, u)

    GT dinv;
    GT::inv(dinv, e2);
    GT::mul(e1, e1, dinv);          // e(v,k1) * e(k0,u)^-1 = e(v,k1) / e(k0,u)
    return e1;                      // this is d_i
}

// BTBF.Combine(pkc=⊥, j, c, J, {d_i}_{i∈J})
BTBF::SymKey BTBF::combine(const std::vector<int>& J, const std::vector<GT>& shares){
    if (J.size() != shares.size()) {
        throw std::invalid_argument("J and shares must have same length");
    }
    if (J.empty()) {
        throw std::invalid_argument("J must be non-empty");
    }

    GT W;
    W.setOne(); // identity in GT

    for (size_t idx = 0; idx < J.size(); idx++) {
        int i = J[idx];                // 1-based index of party
        Fr lambda_i = lagrangeCoeffAtZero(J, i);
        GT term;
        GT::pow(term, shares[idx], lambda_i);
        GT::mul(W, W, term);
    }
    return H2(W);
}

} // namespace btbf
