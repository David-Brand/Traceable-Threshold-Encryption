#include "fingerprinting_codes.hpp"

#include <stdexcept>

namespace fingerprinting {

// ---------- RNG helper ----------

std::mt19937_64& FingerprintingCode::rng() {
    static std::random_device rd;
    static std::mt19937_64 engine(rd());
    return engine;
}

// global RNG for use outside FingerprintingCode (public helper)
std::mt19937_64& global_rng() {
    static std::random_device rd2;
    static std::mt19937_64 engine2(rd2());
    return engine2;
}

// ---------- FingerprintingCode implementation ----------

FingerprintingCode::FingerprintingCode(std::size_t n, double e, std::size_t d_override)
    : n_(n), e_(e)
{
    if (n_ < 2) {
        throw std::invalid_argument("n must be at least 2");
    }
    if (e_ <= 0.0 || e_ >= 1.0) {
        throw std::invalid_argument("e must be in (0,1)");
    }

    if (d_override == 0) {
        // d = ceil(2 * n^2 * log((2n)/e))
        double val = 2.0 * static_cast<double>(n_) * static_cast<double>(n_) * std::log((2.0 * static_cast<double>(n_)) / e_);
        d_ = static_cast<std::size_t>(std::ceil(val));
    } else {
        d_ = d_override;
    }

    l_ = d_ * (n_ - 1);

    // create random permutation
    permutation_.resize(l_);
    std::iota(permutation_.begin(), permutation_.end(), 0);
    std::shuffle(permutation_.begin(), permutation_.end(), rng());

}

void FingerprintingCode::getCodeword(std::size_t user, Bitset& out) const {
    if (out.size() != l_) {
        out.resize(l_, false);
    }
    // fast clear
    std::fill(out.begin(), out.end(), false);

    // same logic as original: set range from user*d_ to l_ (preserve semantics)
    for (std::size_t j = user * d_; j < l_; ++j) {
        out[permutation_[j]] = true;
    }
}

void FingerprintingCode::getCodewordToSlice(std::size_t user, Bitset& out, std::size_t offset) const {

    // same logic as original: set range from user*d_ to l_ (preserve semantics)
    for (std::size_t j = user * d_; j < l_; ++j) {
        out[offset + permutation_[j]] = true;
    }
}

std::size_t FingerprintingCode::weight(const Bitset& x, std::size_t start, std::size_t end) const
{
    if (end > l_) end = l_;
    std::size_t tmp = 0;
    std::size_t idx;
    // local refs for speed
    const auto &perm = permutation_;
    const std::size_t xsz = x.size();
    for (std::size_t i = start; i < end; ++i) {
        idx = perm[i];
        if (idx < xsz && x[idx]) {
            ++tmp;
        }
    }
    return tmp;
}

std::vector<std::size_t> FingerprintingCode::trace(const Bitset& x) const {
    std::vector<std::size_t> guilty;

    if (x.size() < l_) {
        throw std::invalid_argument("trace: input bitset has wrong length");
    }

    // First user
    if (weight(x, 0, d_) > 0) {
        guilty.push_back(0);
    }

    // Middle users: range(1, n-1)
    for (std::size_t i = 1; i + 1 < n_; ++i) {
        std::size_t k = weight(x, (i - 1) * d_, (i + 1) * d_);
        double kd2 = static_cast<double>(k) / 2.0;
        double val = kd2 - std::sqrt(kd2 * std::log((2.0 * static_cast<double>(n_)) / e_));
        std::size_t w = weight(x, (i - 1) * d_, i * d_);
        if (static_cast<double>(w) < val) {
            guilty.push_back(i);
        }
    }

    // Last user
    if (weight(x, l_ - d_, l_) < d_) {
        guilty.push_back(n_ - 1);
    }

    return guilty;
}

Bitset FingerprintingCode::collude(std::size_t i, std::size_t j) const {
    if (i >= n_ || j >= n_) {
        throw std::out_of_range("collude: user index out of range");
    }

    Bitset c(l_, false);
    Bitset a;
    getCodeword(i, a);
    Bitset b;
    getCodeword(j, b);

    for (std::size_t k = 0; k < l_; ++k) {
        c[k] = (a[k] && b[k]);
    }
    return c;
}

// ---------- LogLengthCodes implementation ----------

LogLengthCodes::LogLengthCodes(std::size_t N, std::size_t c, double e)
    : N_(N), c_(c), e_(e)
{
    if (c_ == 0) {
        throw std::invalid_argument("c must be positive");
    }
    if (N_ == 0) {
        throw std::invalid_argument("N must be positive");
    }
    if (e_ <= 0.0 || e_ >= 1.0) {
        throw std::invalid_argument("e must be in (0,1)");
    }

    n_ = 2 * c_;
    std::cout << "n: " << n_ << '\n';

    // L = ceil(2 * c * log((2 * N) / e))
    
    double valL = 2.0 * static_cast<double>(c_) * std::log((2.0 * static_cast<double>(N_)) / e_);
    L_ = static_cast<std::size_t>(std::ceil(valL));
    std::cout << "L: " << L_ << '\n';

    // d = ceil(2 * n^2 * log((4 * n * L) / e))
    
    double valD = 2.0 * static_cast<double>(n_) * static_cast<double>(n_) * std::log((4.0 * static_cast<double>(n_) * static_cast<double>(L_)) / e_);
    d_ = static_cast<std::size_t>(std::ceil(valD));
    std::cout << "d: " << d_ << '\n';

    total_length_ = L_ * d_ * (n_-1);
    std::cout << "length: " << total_length_ << '\n';

    block_length_ = d_ * (n_ - 1);
    std::cout << "block length: " << block_length_ << '\n';

    std::cout << "here\n";

    // Construct L FingerprintingCodes with fixed d_
    components_.reserve(L_);
    for (std::size_t i = 0; i < L_; ++i) {
        components_.emplace_back(n_, e_, d_);
    }

    std::cout << "constructed " << L_ << " component codes\n";

    hiddenCode_ = createHiddenCode();
}

std::vector<std::vector<std::size_t>> LogLengthCodes::createHiddenCode() {
    std::vector<std::vector<std::size_t>> h(N_, std::vector<std::size_t>(L_));
    std::uniform_int_distribution<std::size_t> dist(0, n_ - 1);

    auto& engine = global_rng();
    for (std::size_t i = 0; i < N_; ++i) {
        for (std::size_t j = 0; j < L_; ++j) {
            h[i][j] = dist(engine);
        }
    }
    return h;
}

std::vector<Bitset> LogLengthCodes::getLLCodeword(std::size_t user) const {
    std::vector<Bitset> tmp;
    Bitset w(total_length_, false);
    std::size_t offset;

    for (std::size_t i = 0; i < L_; ++i) {
        offset = i * block_length_;
        components_[i].getCodeword(hiddenCode_[user][i], w);
        tmp.push_back(w);
    }
    return tmp;
}

std::size_t LogLengthCodes::trace(const std::vector<Bitset> x) const {

    if (x.size() != L_ || x[0].size() != block_length_) {
        throw std::invalid_argument("LogLengthCodes::trace: bitset length mismatch");
    }

    // Split x into L blocks and trace each with corresponding code
    std::vector<std::size_t> y;
    y.reserve(L_);

    for (std::size_t i = 0; i < L_; ++i) {
        auto guilty = components_[i].trace(x[i]);
        y.push_back(guilty.empty() ? 0 : guilty[0]);
    }

    // Now find index i with maximum matches between hiddenCode_[i] and y
    std::size_t bestIndex = 0;
    std::size_t bestMatches = 0;

    for (std::size_t i = 0; i < N_; ++i) {
        std::size_t count = 0;
        for (std::size_t j = 0; j < L_; ++j) {
            if (hiddenCode_[i][j] == y[j]) {
                ++count;
            }
        }
        if (count > bestMatches) {
            bestMatches = count;
            bestIndex = i;
        }
    }
    return bestIndex;
}

std::vector<Bitset> LogLengthCodes::collude(const std::vector<std::size_t>& coalition) const {

    std::vector<Bitset> word;

    for(std::size_t i = 0; i < L_; ++i){
        word.push_back(Bitset(block_length_, false));
    }

    std::vector<Bitset> vec;

    for(auto idx : coalition){
        if (idx >= N_) {
            throw std::out_of_range("LogLengthCodes::collude: coalition index out of range");
        }

        vec = getLLCodeword(idx);
        for(std::size_t j = 0; j < L_; ++j){
            for (std::size_t k = 0; k < block_length_; ++k) {
                if(word[j][k] || !vec[j][k]) continue; // already set to one or tmp is zero
                word[j][k] = vec[j][k];
            }
        }
    }

    return word;
}

} // namespace fingerprinting
