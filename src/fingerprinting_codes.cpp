// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: MIT

#include "fingerprinting_codes.hpp"
#include <string>
#include <algorithm>
#include <cmath>


namespace fingerprinting {


// ---------- FingerprintingCode implementation ----------

// Constructs a fingerprinting code with parameters n, e, delta.
// Computes block and word lengths, optionally using an override, and stores RNG/permutation.
FingerprintingCode::FingerprintingCode(std::size_t n, double e, float delta, std::mt19937_64* engine, std::vector<std::size_t>* permutation, std::size_t d_override) : num_fingerprints_(n), tracing_error_(e), engine_(engine), permutation_(permutation) {
    // Validate inputs
    if (num_fingerprints_ < 2) throw std::invalid_argument("n must be at least 2");
    if (tracing_error_ <= 0.0 || tracing_error_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");

    // Compute block length (or use override)
    if (d_override == 0) {
        double val = ((4 * std::pow(num_fingerprints_, 2))/std::pow(1-delta, 2)) * std::log((2.0 * num_fingerprints_) / tracing_error_);
        block_length_ = static_cast<std::size_t>(std::ceil(val));
    } else {
        block_length_ = d_override;
    }
    // Derive total word length
    word_length_ = block_length_ * (num_fingerprints_ - 1);
}

// Produces the codeword for a given user by setting bits via a permutation.
// The output bitset is resized and populated using the stored permutation.
void FingerprintingCode::getCodeword(std::size_t user, PackedBitset& out) const {
    // Prepare output bitset
    out.resize(word_length_);
    out.clear();
    // Set bits for the user's block via permutation
    for (std::size_t j = user * block_length_; j < word_length_; ++j) {
        out.set((*permutation_)[j]);
    }
}

// Counts the number of set bits in x between start and end using the permutation.
// Returns the weight of that segment.
std::size_t FingerprintingCode::weight(const PackedBitset& x, std::size_t start, std::size_t end) const {
    // Clamp end to word length
    if (end > word_length_) end = word_length_;
    std::size_t tmp = 0;
    // Accumulate number of set bits in the range
    for (std::size_t i = start; i < end; ++i) {
        if (x.get((*permutation_)[i])) ++tmp;
    }
    return tmp;
}

// Traces guilty users by evaluating segment weights against thresholds.
// Returns indices suspected of participating in the collusion.
std::vector<std::size_t> FingerprintingCode::trace(const PackedBitset& x, const PackedBitset& x_unreadable) const {
    std::vector<std::size_t> guilty;
    // Validate bitset length
    if (x.size() != word_length_) throw std::invalid_argument("trace: input bitset has wrong length");

    // Check first and last users via boundary blocks
    if ((weight(x, 0, block_length_) > 0) || (weight(x_unreadable, 0, block_length_) > 0)) guilty.push_back(0);

    // Sweep middle users and compare weights against thresholds
    for (std::size_t i = 1; i + 1 < num_fingerprints_; ++i) {
        std::size_t k = weight(x, (i - 1) * block_length_, (i + 1) * block_length_);
        double kd2 = static_cast<double>(k) / 2.0;
        double val = kd2 - std::sqrt(kd2 * std::log((2.0 * num_fingerprints_) / tracing_error_));
        std::size_t w = weight(x, (i - 1) * block_length_, i * block_length_);
        if (static_cast<double>(w) < val){ 
            guilty.push_back(i);
            continue;
        }

        std::size_t k_unreadable = weight(x_unreadable, (i - 1) * block_length_, (i + 1) * block_length_);
        double kd2_unreadable = static_cast<double>(k_unreadable) / 2.0;
        double val_unreadable = kd2_unreadable - std::sqrt(kd2_unreadable * std::log((2.0 * num_fingerprints_) / tracing_error_));
        std::size_t w_unreadable = weight(x_unreadable, (i - 1) * block_length_, i * block_length_);
        if (static_cast<double>(w_unreadable) < val_unreadable){ 
            guilty.push_back(i);
        }
    }

    // Final boundary checks
    if ((weight(x, word_length_ - block_length_, word_length_) < block_length_) || (weight(x_unreadable, word_length_ - block_length_, word_length_) > 0)) guilty.push_back(num_fingerprints_ - 1);
    
    return guilty;
}

// OR-colludes a coalition’s codewords to form an illegal word.
// Returns the combined bitset.
PackedBitset FingerprintingCode::collude(const std::vector<std::size_t>& coalition) const {
    // Initialize accumulator bitset
    PackedBitset a(word_length_);

    // Combine all colluders’ codewords by OR
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(coalition.size()); ++i){
        PackedBitset cw(word_length_);
        getCodeword(coalition[i], cw);
        a.orInPlace(cw);
    }
    return a;
}

// ---------- LogLengthCodes implementation ----------

// Constructs log-length codes with parameters for user count, coalition threshold, tracing error, and delta.
// Initializes RNG and computes necessary lengths and component counts, adjusting for memory limits.
LogLengthCodes::LogLengthCodes(std::size_t N, std::size_t c, double e, float delta) : NumFingerprints_(N), coalition_threshold_(c), tracing_error_(e) {
    // Validate inputs
    if (coalition_threshold_ == 0) throw std::invalid_argument("c must be positive");
    if (NumFingerprints_ == 0) throw std::invalid_argument("N must be positive");
    if (tracing_error_ <= 0.0 || tracing_error_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");

    std::random_device rd;
    engine_ = std::mt19937_64(rd());

    //initialize parameters so e holds
    word_length_ = 2 * coalition_threshold_;
    double valL = 2.0 * coalition_threshold_ * std::log((2.0 * NumFingerprints_) / tracing_error_);
    component_count_ = static_cast<std::size_t>(std::ceil(valL));
    double valD = ((4 * std::pow(word_length_, 2))/std::pow(1-delta, 2)) * std::log((4.0 * word_length_ * component_count_) / tracing_error_);
    block_length_ = static_cast<std::size_t>(std::ceil(valD));
    total_length_ = component_count_ * block_length_ * (word_length_-1);
    component_length_ = block_length_ * (word_length_ - 1);

    double tmpd = ((4 * std::pow(NumFingerprints_, 2))/std::pow(1-delta, 2)) * std::log((2.0 * NumFingerprints_) / tracing_error_);
    if (total_length_ > static_cast<std::size_t>(std::ceil(tmpd))*(NumFingerprints_-1)){
        component_count_ = 1;
        block_length_ = static_cast<std::size_t>(std::ceil(tmpd));
        word_length_ = NumFingerprints_;
        total_length_ = component_count_ * block_length_ * (word_length_-1);
        component_length_ = total_length_;
    }

    // error if memory usage would be too high
    uint64_t totalBits = uint64_t(block_length_) * uint64_t(word_length_ - 1);
    double bytes = double(totalBits) / 8.0;
    const double GiB = 1024.0 * 1024.0 * 1024.0;
    if (bytes > 14.0 * GiB)
        throw std::runtime_error("LogLengthCodes: allocation too large (" + std::to_string(bytes / GiB) + " GiB). Aborting.");


    //compute permutations, bound to 10 repeating instances to save memory
    std::size_t pCount = component_count_ > 10 ? 10 : component_count_;
    permutations_.resize(pCount);
    for(auto& perm : permutations_){
        perm.resize(component_length_);
        std::iota(perm.begin(), perm.end(), 0);
        std::shuffle(perm.begin(), perm.end(), engine_);
    }

    //component initialization
    components_.reserve(component_count_);
    for (std::size_t i = 0; i < component_count_; ++i)
        components_.emplace_back(word_length_, tracing_error_, delta, &engine_, &permutations_[i%10], block_length_);
    hiddenCode_ = createHiddenCode();
}

// Creates the hidden code matrix for all users.
// Fills with random indices for each user's codeword components.
std::vector<std::vector<std::size_t>> LogLengthCodes::createHiddenCode() {
    std::vector<std::vector<std::size_t>> h(NumFingerprints_, std::vector<std::size_t>(component_count_));
    std::uniform_int_distribution<std::size_t> dist(0, word_length_ - 1);

    //fill hidden code with random numbers from 0 to n-1
    for (std::size_t i = 0; i < NumFingerprints_; ++i)
        for (std::size_t j = 0; j < component_count_; ++j)
            h[i][j] = dist(engine_);
    return h;
}

// Returns the set of component bitsets representing one user’s log-length codeword.
// Fills out with the codeword for the hidden index in each component.
void LogLengthCodes::getLLCodeword(std::size_t user, std::vector<PackedBitset>& out) const {
    // Resize output to component count
    out.resize(component_count_);
    // Fetch each component codeword for the hidden selection
    for (std::size_t i = 0; i < component_count_; ++i)
        components_[i].getCodeword(hiddenCode_[user][i], out[i]);
}

// Colludes a coalition’s codewords by computing the OR across components.
// Returns the combined codewords for the coalition.
std::vector<PackedBitset> LogLengthCodes::collude(const std::vector<std::size_t>& coalition) const {
    std::vector<PackedBitset> out(component_count_, PackedBitset(component_length_));

    // Parallelize across components
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(component_count_); ++i) {

        PackedBitset cw(component_length_);  // thread-local buffer

        for (auto idx : coalition) {
            //codeword for user idx in component i
            components_[i].getCodeword(hiddenCode_[idx][i], cw);

            // Block-wise OR
            out[i].orInPlace(cw);
        }
    }
    return out;
}


// Traces the source of a codeword by comparing with all hidden codewords.
// Returns the index of the user whose codeword matches best.
std::size_t LogLengthCodes::trace(const std::vector<PackedBitset>& x, const std::vector<PackedBitset>& x_unreadable) const {
    if (x.size() != component_count_ || x[0].size() != component_length_)
        throw std::invalid_argument("LogLengthCodes::trace: bitset length mismatch");

    std::vector<std::size_t> y(component_count_);

    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(component_count_); ++i) {
        auto guilty = components_[i].trace(x[static_cast<std::size_t>(i)], x_unreadable[static_cast<std::size_t>(i)]);
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
        for (long long i = 0; i < static_cast<long long>(NumFingerprints_); ++i) {
            std::size_t user = static_cast<std::size_t>(i);
            const auto& row = hiddenCode_[user];

            std::size_t count = 0;          //thread-local match count
            for (std::size_t j = 0; j < component_count_; ++j) {
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

//------------ Tardos Codes implementation -----------

// Initializes Tardos codes with coalition threshold and error.
// Precomputes probabilities, codebook, and thresholds used for tracing.
TardosCodes::TardosCodes(std::size_t c, double e, std::size_t n) : coalition_threshold_(c), tracing_error_(e), num_fingerprints_(n) {
    // Validate inputs and precompute parameters
    if (coalition_threshold_ == 0) throw std::invalid_argument("c must be positive");
    if (tracing_error_ <= 0.0 || tracing_error_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");
    probabilities_ = calculateProbabilities(coalition_threshold_, tracing_error_, num_fingerprints_);
    auto tmp = generateCodeBook(probabilities_, num_fingerprints_);
    codeBook_ = tmp.first;
    U_ = tmp.second;
    // Compute code length and thresholds
    k_ = std::log(static_cast<double>(std::max<std::size_t>(static_cast<std::size_t>(1), num_fingerprints_)) / tracing_error_);
    word_length_ = static_cast<std::size_t>(std::ceil(100 * std::pow(coalition_threshold_, 2) * k_));
    accusation_threshold_ = 20 * coalition_threshold_ * k_;
}

std::vector<double> TardosCodes::calculateProbabilities(std::size_t cs, double err, std::size_t new_words){
    std::size_t m = static_cast<std::size_t>(std::ceil(100 * std::pow(cs, 2) * std::log(static_cast<double>(std::max<std::size_t>(static_cast<std::size_t>(1), new_words)) / err)));
    double t = 1.0 / (300.0 * static_cast<double>(cs));
    double tPrime = std::asin(std::sqrt(t));
    double pi = 4.0 * std::atan(1.0);
    if(tPrime >= pi/4.0) tPrime = pi/4.00001;
    if(tPrime <= 0.0) tPrime = 0.00000000001;
    
    // bounds for r_i : [t', π/2 - t']
    double left = tPrime;
    double right = (pi / 2.0) - tPrime;

    std::mt19937_64 random(std::random_device{}());
    std::uniform_real_distribution<double> dist(left, right);

    std::vector<double> probabilities;
    probabilities.reserve(m);
    for(std::size_t i = 0; i < m; ++i){
        double r = dist(random);
        probabilities.push_back(std::pow(std::sin(r), 2.0));
    }
    return probabilities;
}

std::pair<PackedBitset, std::vector<double>> TardosCodes::writeCodeWord(){
    // Generate and append a fresh codeword
    auto tmp = generateCodeWord(probabilities_);
    codeBook_.push_back(tmp.first);
    U_.push_back(tmp.second);
    // Increase user count and return
    num_fingerprints_++;
    return tmp;
}

// Generates a single Tardos codeword from probabilities.
// Also returns the corresponding U-weights for tracing.
std::pair<PackedBitset, std::vector<double>> TardosCodes::generateCodeWord(const std::vector<double>& probabilities){
    // Sample bits according to probabilities
    std::mt19937_64 random(std::random_device{}());
    std::uniform_real_distribution<double> dist(0, 1);

    // Build codeword and weights
    PackedBitset codeword(probabilities.size());
    std::vector<double> U(probabilities.size());
    for(std::size_t i = 0; i < probabilities.size(); ++i){
        double rnd = dist(random);
        if(rnd < probabilities[i]){
            codeword.set(i, true);
            U[i] = std::sqrt((1.0 - probabilities[i]) / probabilities[i]);
        } else {
            codeword.set(i, false);
            U[i] = -std::sqrt(probabilities[i] / (1.0 - probabilities[i]));
        }
    }
    return {codeword, U};
}

// Appends a batch of new codewords to the codebook.
// Returns the generated codewords and U-weights.
std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> TardosCodes::writeCodeBook(std::size_t new_words){
    // Generate a batch of codewords
    auto tmp = generateCodeBook(probabilities_, new_words);
    // Append to existing codebook/state
    codeBook_.insert(codeBook_.end(), tmp.first.begin(), tmp.first.end());
    U_.insert(U_.end(), tmp.second.begin(), tmp.second.end());
    num_fingerprints_ += new_words;
    return tmp;
}

// Generates a batch of Tardos codewords from probabilities.
// Returns all codewords and their U-weights.
std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> TardosCodes::generateCodeBook(const std::vector<double>& probabilities, std::size_t new_words){
    // Prepare containers
    std::vector<PackedBitset> codeBook;
    std::vector<std::vector<double>> Ulist;
    codeBook.reserve(new_words);
    Ulist.reserve(new_words);
    // Generate each codeword independently
    for(std::size_t i = 0; i < new_words; ++i){
        auto tmp = generateCodeWord(probabilities);
        codeBook.push_back(tmp.first);
        Ulist.push_back(tmp.second);
    }
    return {codeBook, Ulist};
}

// Traces suspects using the precomputed U-weights and accusation threshold.
// Delegates to the generic trace overload.
std::vector<std::size_t> TardosCodes::trace(PackedBitset& y){
    // Use stored U and threshold for tracing
    return trace(U_, y, accusation_threshold_);
}

// Scores each user with U-weights against the observed y and threshold z.
// Returns indices of accused users.
std::vector<std::size_t> TardosCodes::trace(const std::vector<std::vector<double>>& U, PackedBitset& y, double z){
    std::vector<std::size_t> accused;
    // Compute matching score per user
    for(std::size_t j = 0; j < U.size(); ++j){
        double sum = 0.0;
        for(std::size_t i = 0; i < y.size(); ++i){
            if(y.get(i)){
                sum += U[j][i];
            }
        }
        // Accuse if score exceeds threshold
        if(sum > z){
            accused.push_back(j);
        }
    }
    return accused;
}

// Colludes codewords by OR to form a suspect word.
// Returns the combined bitset of the coalition.
PackedBitset TardosCodes::collude(const std::vector<std::size_t>& coalition) const{
    // Aggregate all colluders’ codewords
    PackedBitset out(word_length_);
    for(auto idx : coalition){
        out.orInPlace(codeBook_[idx]);
    }
    return out;
}

} // namespace fingerprinting
