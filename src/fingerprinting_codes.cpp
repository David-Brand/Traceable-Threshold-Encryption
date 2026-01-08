#include "fingerprinting_codes.hpp"
#include <string>
#include <algorithm>
#include <cmath>


namespace fingerprinting {


// ---------- FingerprintingCode implementation ----------

FingerprintingCode::FingerprintingCode(std::size_t n, double e, float delta, std::mt19937_64* engine, std::vector<std::size_t>* permutation, std::size_t d_override) : num_fingerprints_(n), tracing_error_(e), engine_(engine), permutation_(permutation) {
    if (num_fingerprints_ < 2) throw std::invalid_argument("n must be at least 2");
    if (tracing_error_ <= 0.0 || tracing_error_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");

    if (d_override == 0) {
        double val = ((4 * std::pow(num_fingerprints_, 2))/std::pow(1-delta, 2)) * std::log((2.0 * num_fingerprints_) / tracing_error_);
        block_length_ = static_cast<std::size_t>(std::ceil(val));
    } else {
        block_length_ = d_override;
    }
    word_length_ = block_length_ * (num_fingerprints_ - 1);
}

void FingerprintingCode::getCodeword(std::size_t user, PackedBitset& out) const {
    out.resize(word_length_);
    out.clear();
    for (std::size_t j = user * block_length_; j < word_length_; ++j) {
        out.set((*permutation_)[j]);
    }
}

std::size_t FingerprintingCode::weight(const PackedBitset& x, std::size_t start, std::size_t end) const {
    if (end > word_length_) end = word_length_;
    std::size_t tmp = 0;
    for (std::size_t i = start; i < end; ++i) {
        if (x.get((*permutation_)[i])) ++tmp;
    }
    return tmp;
}

std::vector<std::size_t> FingerprintingCode::trace(const PackedBitset& x, const PackedBitset& x_unreadable) const {
    std::vector<std::size_t> guilty;
    if (x.size() != word_length_) throw std::invalid_argument("trace: input bitset has wrong length");

    if ((weight(x, 0, block_length_) > 0) || (weight(x_unreadable, 0, block_length_) > 0)) guilty.push_back(0);

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

    if ((weight(x, word_length_ - block_length_, word_length_) < block_length_) || (weight(x_unreadable, word_length_ - block_length_, word_length_) > 0)) guilty.push_back(num_fingerprints_ - 1);
    
    return guilty;
}

PackedBitset FingerprintingCode::collude(const std::vector<std::size_t>& coalition) const {
    PackedBitset a(word_length_);

    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < static_cast<long long>(coalition.size()); ++i){
        PackedBitset cw(word_length_);
        getCodeword(coalition[i], cw);
        a.orInPlace(cw);
    }
    return a;
}

// ---------- LogLengthCodes implementation ----------

LogLengthCodes::LogLengthCodes(std::size_t N, std::size_t c, double e, float delta) : NumFingerprints_(N), coalition_threshold_(c), tracing_error_(e) {
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

std::vector<std::vector<std::size_t>> LogLengthCodes::createHiddenCode() {
    std::vector<std::vector<std::size_t>> h(NumFingerprints_, std::vector<std::size_t>(component_count_));
    std::uniform_int_distribution<std::size_t> dist(0, word_length_ - 1);

    //fill hidden code with random numbers from 0 to n-1
    for (std::size_t i = 0; i < NumFingerprints_; ++i)
        for (std::size_t j = 0; j < component_count_; ++j)
            h[i][j] = dist(engine_);
    return h;
}

void LogLengthCodes::getLLCodeword(std::size_t user, std::vector<PackedBitset>& out) const {
    out.resize(component_count_);
    for (std::size_t i = 0; i < component_count_; ++i)
        components_[i].getCodeword(hiddenCode_[user][i], out[i]);
}

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

TardosCodes::TardosCodes(std::size_t c, double e, std::size_t n) : coalition_threshold_(c), tracing_error_(e), num_fingerprints_(n) {
    if (coalition_threshold_ == 0) throw std::invalid_argument("c must be positive");
    if (tracing_error_ <= 0.0 || tracing_error_ >= 1.0) throw std::invalid_argument("e must be in (0,1)");
    probabilities_ = calculateProbabilities(coalition_threshold_, tracing_error_, num_fingerprints_);
    auto tmp = generateCodeBook(probabilities_, num_fingerprints_);
    codeBook_ = tmp.first;
    U_ = tmp.second;
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
    auto tmp = generateCodeWord(probabilities_);
    codeBook_.push_back(tmp.first);
    U_.push_back(tmp.second);
    num_fingerprints_++;
    return tmp;
}

std::pair<PackedBitset, std::vector<double>> TardosCodes::generateCodeWord(const std::vector<double>& probabilities){
    std::mt19937_64 random(std::random_device{}());
    std::uniform_real_distribution<double> dist(0, 1);

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

std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> TardosCodes::writeCodeBook(std::size_t new_words){
    auto tmp = generateCodeBook(probabilities_, new_words);
    codeBook_.insert(codeBook_.end(), tmp.first.begin(), tmp.first.end());
    U_.insert(U_.end(), tmp.second.begin(), tmp.second.end());
    num_fingerprints_ += new_words;
    return tmp;
}

std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> TardosCodes::generateCodeBook(const std::vector<double>& probabilities, std::size_t new_words){
    std::vector<PackedBitset> codeBook;
    std::vector<std::vector<double>> Ulist;
    codeBook.reserve(new_words);
    Ulist.reserve(new_words);
    for(std::size_t i = 0; i < new_words; ++i){
        auto tmp = generateCodeWord(probabilities);
        codeBook.push_back(tmp.first);
        Ulist.push_back(tmp.second);
    }
    return {codeBook, Ulist};
}

std::vector<std::size_t> TardosCodes::trace(PackedBitset& y){
    return trace(U_, y, accusation_threshold_);
}

std::vector<std::size_t> TardosCodes::trace(const std::vector<std::vector<double>>& U, PackedBitset& y, double z){
    std::vector<std::size_t> accused;
    for(std::size_t j = 0; j < U.size(); ++j){
        double sum = 0.0;
        for(std::size_t i = 0; i < y.size(); ++i){
            if(y.get(i)){
                sum += U[j][i];
            }
        }
        if(sum > z){
            accused.push_back(j);
        }
    }
    return accused;
}

PackedBitset TardosCodes::collude(const std::vector<std::size_t>& coalition) const{
    PackedBitset out(word_length_);
    for(auto idx : coalition){
        out.orInPlace(codeBook_[idx]);
    }
    return out;
}

} // namespace fingerprinting
