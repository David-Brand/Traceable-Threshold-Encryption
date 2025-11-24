#include "fingerprinting_codes.hpp"
#include <string>

namespace fingerprinting {


// ---------- FingerprintingCode implementation ----------

FingerprintingCode::FingerprintingCode(std::size_t n, double e, std::mt19937_64& engine, std::vector<std::size_t>& permutation, std::size_t d_override)
    : n_(n), e_(e), engine_(engine), permutation_(permutation)
{
    if (n_ < 2) throw std::invalid_argument("n must be at least 2");
    if (e_ <= 0.0 || e_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");

    if (d_override == 0) {
        double val = 2.0 * n_ * n_ * std::log((2.0 * n_) / e_);
        d_ = static_cast<std::size_t>(std::ceil(val));
    } else {
        d_ = d_override;
    }
    l_ = d_ * (n_ - 1);
}

void FingerprintingCode::getCodeword(std::size_t user, PackedBitset& out) const {
    out.resize(l_);
    out.clear();
    for (std::size_t j = user * d_; j < l_; ++j) {
        out.set(permutation_[j]);
    }
}

std::size_t FingerprintingCode::weight(const PackedBitset& x, std::size_t start, std::size_t end) const {
    if (end > l_) end = l_;
    std::size_t tmp = 0;
    for (std::size_t i = start; i < end; ++i) {
        if (x.get(permutation_[i])) ++tmp;
    }
    return tmp;
}

std::vector<std::size_t> FingerprintingCode::trace(const PackedBitset& x) const {
    std::vector<std::size_t> guilty;
    if (x.size() != l_) throw std::invalid_argument("trace: input bitset has wrong length");

    if (weight(x, 0, d_) > 0) guilty.push_back(0);

    for (std::size_t i = 1; i + 1 < n_; ++i) {
        std::size_t k = weight(x, (i - 1) * d_, (i + 1) * d_);
        double kd2 = static_cast<double>(k) / 2.0;
        double val = kd2 - std::sqrt(kd2 * std::log((2.0 * n_) / e_));
        std::size_t w = weight(x, (i - 1) * d_, i * d_);
        if (static_cast<double>(w) < val) guilty.push_back(i);
    }

    if (weight(x, l_ - d_, l_) < d_) guilty.push_back(n_ - 1);
    
    return guilty;
}

void FingerprintingCode::collude(std::size_t i, std::size_t j, PackedBitset& out) const {
    PackedBitset a(l_), b(l_);
    getCodeword(i, a);
    getCodeword(j, b);
    out.resize(l_);
    out.clear();
    for (std::size_t k = 0; k < l_; ++k)
        out.set(k, a.get(k) || b.get(k));
}

// ---------- LogLengthCodes implementation ----------

LogLengthCodes::LogLengthCodes(std::size_t N, std::size_t c, double e)
    : N_(N), c_(c), e_(e)
{
    if (c_ == 0) throw std::invalid_argument("c must be positive");
    if (N_ == 0) throw std::invalid_argument("N must be positive");
    if (e_ <= 0.0 || e_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");

    std::random_device rd;
    engine_ = std::mt19937_64(rd());

    //initialize parameters so e holds
    n_ = 2 * c_;
    double valL = 2.0 * c_ * std::log((2.0 * N_) / e_);
    L_ = static_cast<std::size_t>(std::ceil(valL));
    double valD = 2.0 * n_ * n_ * std::log((4.0 * n_ * L_) / e_);
    d_ = static_cast<std::size_t>(std::ceil(valD));
    total_length_ = L_ * d_ * (n_-1);
    block_length_ = d_ * (n_ - 1);

    // error if memory usage would be too high
    uint64_t totalBits = uint64_t(L_) * uint64_t(d_) * uint64_t(n_ - 1);
    double bytes = double(totalBits) / 8.0;
    const double GiB = 1024.0 * 1024.0 * 1024.0;
    if (bytes > 8.0 * GiB)
        throw std::runtime_error("LogLengthCodes: allocation too large (" + std::to_string(bytes / GiB) + " GiB). Aborting.");

    
    //compute permutations, bound to 100 repeating instances to save memory
    permutations_.resize(100);
    for(auto& permutation_ : permutations_) {
        permutation_.resize(block_length_);
        std::iota(permutation_.begin(), permutation_.end(), 0);
        std::shuffle(permutation_.begin(), permutation_.end(), engine_);
    }

    //component initialization
    components_.reserve(L_);
    for (std::size_t i = 0; i < L_; ++i)
        components_.emplace_back(n_, e_, engine_, permutations_[i%100], d_);
    hiddenCode_ = createHiddenCode();
}

std::vector<std::vector<std::size_t>> LogLengthCodes::createHiddenCode() {
    std::vector<std::vector<std::size_t>> h(N_, std::vector<std::size_t>(L_));
    std::uniform_int_distribution<std::size_t> dist(0, n_ - 1);

    //fill hidden code with random numbers from 0 to n-1
    for (std::size_t i = 0; i < N_; ++i)
        for (std::size_t j = 0; j < L_; ++j)
            h[i][j] = dist(engine_);
    return h;
}

void LogLengthCodes::getLLCodeword(std::size_t user, std::vector<PackedBitset>& out) const {
    out.resize(L_);
    for (std::size_t i = 0; i < L_; ++i)
        components_[i].getCodeword(hiddenCode_[user][i], out[i]);
}

std::vector<PackedBitset> LogLengthCodes::collude(
    const std::vector<std::size_t>& coalition) const
{
    std::vector<PackedBitset> out(L_, PackedBitset(block_length_));

    // Parallelize across components, static schedule because each iteration has about same workload
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(L_); ++i) {

        PackedBitset cw(block_length_);  // thread-local buffer

        for (auto idx : coalition) {
            //codeword for user idx in component i
            components_[i].getCodeword(hiddenCode_[idx][i], cw);

            // Block-wise OR
            out[i].orInPlace(cw);
        }
    }
    return out;
}


std::size_t LogLengthCodes::trace(const std::vector<PackedBitset>& x) const {
    if (x.size() != L_ || x[0].size() != block_length_)
        throw std::invalid_argument("LogLengthCodes::trace: bitset length mismatch");

    std::vector<std::size_t> y(L_);

    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(L_); ++i) {
        auto guilty = components_[i].trace(x[static_cast<std::size_t>(i)]);
        y[static_cast<std::size_t>(i)] = guilty.empty() ? 0 : guilty[0];
    }
    std::size_t bestIndex   = 0;
    std::size_t bestMatches = 0;

    // find word with most matches with y
    #pragma omp parallel
    {
        std::size_t localBestIndex   = 0;
        std::size_t localBestMatches = 0;

        #pragma omp for nowait schedule(static)
        for (long long i = 0; i < static_cast<long long>(N_); ++i) {
            std::size_t user = static_cast<std::size_t>(i);
            const auto& row = hiddenCode_[user];

            std::size_t count = 0;          //thread-local match count
            for (std::size_t j = 0; j < L_; ++j) {
                if (row[j] == y[j]) ++count;
            }

            if (count > localBestMatches) {
                localBestMatches = count;
                localBestIndex   = user;
            }
        }
        // Combine thread-local bests into global best
        #pragma omp critical
        {
            if (localBestMatches > bestMatches) {
                bestMatches = localBestMatches;
                bestIndex   = localBestIndex;
            }
        }
    }
    return bestIndex;
}

} // namespace fingerprinting
