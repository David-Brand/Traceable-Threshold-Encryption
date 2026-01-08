#ifndef FINGERPRINTING_CODES_HPP
#define FINGERPRINTING_CODES_HPP

#include <vector>
#include <cstddef>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <intrin.h>

namespace fingerprinting {

// Packed bitset: stores bits in uint64_t blocks
class PackedBitset {
public:
    PackedBitset(std::size_t n = 0) : bits_((n+63)/64, 0), size_(n) {}
    void resize(std::size_t n) {
        bits_.resize((n+63)/64, 0);
        size_ = n;
    }
    void clear() { std::fill(bits_.begin(), bits_.end(), 0); }
    bool get(std::size_t i) const {
        return (bits_[i/64] >> (i%64)) & 1;
    }
    void set(std::size_t i, bool v = true) {
        if (v)
            bits_[i/64] |= (uint64_t(1) << (i%64));
        else
            bits_[i/64] &= ~(uint64_t(1) << (i%64));
    }
    std::size_t size() const { return size_; }

    void orInPlace(const PackedBitset& other) {
        auto& a = bits_;
        const auto& b = other.blocks();
        std::size_t n = a.size();
        for (std::size_t i = 0; i < n; ++i) {
            a[i] |= b[i];
        }
    }


    // direct access to blocks
    const std::vector<uint64_t>& blocks() const { return bits_; }
    std::vector<uint64_t>& blocks() { return bits_; }
private:
    std::vector<uint64_t> bits_;
    std::size_t size_;
};

class FingerprintingCode {
public:
    FingerprintingCode(std::size_t num_fingerprints, double tracing_error, float delta, std::mt19937_64* engine, std::vector<std::size_t>* permutation, std::size_t d_override = 0);

    std::size_t get_num_fingerprints() const { return num_fingerprints_; }
    double get_tracing_error() const { return tracing_error_; }
    std::size_t get_block_length() const { return block_length_; }
    std::size_t get_word_length() const { return word_length_; }
    std::vector<std::size_t> get_permutation() const { return *permutation_;}

    std::size_t getLength() const {return word_length_; }

    // Generate codeword for user
    void getCodeword(std::size_t user, PackedBitset& out) const;

    //returns indices of guilty users
    std::vector<std::size_t> trace(const PackedBitset& x, const PackedBitset& x_unreadable) const;

    PackedBitset collude(const std::vector<std::size_t>& coalition) const;

private:
    std::size_t num_fingerprints_, block_length_, word_length_;
    double tracing_error_;
    std::vector<std::size_t>* permutation_;
    std::mt19937_64* engine_;

    std::size_t weight(const PackedBitset& x, std::size_t start, std::size_t end) const;
};

class LogLengthCodes {
public:
    LogLengthCodes(std::size_t NumFingerprints, std::size_t coalition_threshold, double tracing_error, float delta);

    std::size_t NumFingerprints() const { return NumFingerprints_; }
    std::size_t coalition_threshold() const { return coalition_threshold_; }
    std::size_t component_count() const { return component_count_; }
    std::size_t word_length() const { return word_length_; }
    std::size_t block_length() const { return block_length_; }
    std::size_t getLength() const {return total_length_; }
    std::size_t componentLength() const {return component_length_; }
    std::vector<FingerprintingCode> components() const { return components_;}

    //return index of most likely guilty word
    std::size_t trace(const std::vector<PackedBitset>& x, const std::vector<PackedBitset>& x_unreadable) const;

    std::vector<PackedBitset> collude(const std::vector<std::size_t>& coalition) const;

    // Generate codeword for user
    void getLLCodeword(std::size_t user, std::vector<PackedBitset>& out) const;

private:
    std::size_t NumFingerprints_, coalition_threshold_, component_count_, word_length_, block_length_, total_length_, component_length_;
    double tracing_error_;
    std::vector<FingerprintingCode> components_;
    std::vector<std::vector<std::size_t>> hiddenCode_;
    std::mt19937_64 engine_;
    std::vector<std::vector<std::size_t>> permutations_;

    std::vector<std::vector<std::size_t>> createHiddenCode();
};

class TardosCodes{
    public:
    TardosCodes(std::size_t coalition_threshold, double tracing_error, std::size_t num_fingerprints = 0);

    static std::pair<PackedBitset, std::vector<double>> generateCodeWord(const std::vector<double>& probabilities);
    static std::vector<double> calculateProbabilities(std::size_t cs, double err, std::size_t new_words=0);
    static std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> generateCodeBook(const std::vector<double>& probabilities, std::size_t new_words);
    static std::vector<std::size_t> trace(const std::vector<std::vector<double>>& U, PackedBitset& y, double z);

    std::pair<PackedBitset, std::vector<double>> writeCodeWord();
    std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> writeCodeBook(std::size_t new_words);
    std::vector<std::size_t> trace(PackedBitset& y);
    PackedBitset collude(const std::vector<std::size_t>& coalition) const;
    bool getBit(std::size_t user, std::size_t pos) const {
        return codeBook_[user].get(pos);
    }

    std::size_t get_coalition_threshold() const { return coalition_threshold_; }
    std::size_t get_num_fingerprints() const { return num_fingerprints_; }
    double get_tracing_error() const { return tracing_error_; }
    double k() const { return k_; }
    double get_accusation_threshold() const { return accusation_threshold_; }
    std::size_t getLength() const { return word_length_; }
    std::vector<double> getProbabilities() const { return probabilities_;}
    std::vector<PackedBitset> getCodeBook() const { return codeBook_; }
    std::vector<std::vector<double>> getU() const { return U_; }

    private:
    std::vector<double> probabilities_;
    std::vector<PackedBitset> codeBook_;
    std::vector<std::vector<double>> U_;

    std::size_t word_length_, coalition_threshold_, num_fingerprints_;
    double tracing_error_, k_, accusation_threshold_;
};

} // namespace fingerprinting

#endif // FINGERPRINTING_CODES_HPP
