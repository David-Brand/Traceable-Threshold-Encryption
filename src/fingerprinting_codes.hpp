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
    FingerprintingCode(std::size_t n, double e, float delta, std::mt19937_64& engine, std::vector<std::size_t>& permutation, std::size_t d_override = 0);

    std::size_t n() const { return n_; }
    double e() const { return e_; }
    std::size_t d() const { return d_; }
    std::size_t l() const { return l_; }
    std::vector<std::size_t> permutation() const { return permutation_;}

    // Generate codeword for user
    void getCodeword(std::size_t user, PackedBitset& out) const;

    //returns indices of guilty users
    std::vector<std::size_t> trace(const PackedBitset& x, const PackedBitset& x_unreadable) const;

    //OR of two codewords
    void collude(std::size_t i, std::size_t j, PackedBitset& out) const;

private:
    std::size_t n_, d_, l_;
    double e_;
    std::vector<std::size_t>& permutation_;
    std::mt19937_64& engine_;

    std::size_t weight(const PackedBitset& x, std::size_t start, std::size_t end) const;
};

class LogLengthCodes {
public:
    LogLengthCodes(std::size_t N, std::size_t c, double e, float delta);

    std::size_t N() const { return N_; }
    std::size_t c() const { return c_; }
    std::size_t L() const { return L_; }
    std::size_t n() const { return n_; }
    std::size_t d() const { return d_; }
    std::size_t totalLength() const {return total_length_; }
    std::size_t blockLength() const {return block_length_; }
    std::vector<FingerprintingCode> components() const { return components_;}

    //return index of most likely guilty word
    std::size_t trace(const std::vector<PackedBitset>& x, const std::vector<PackedBitset>& x_unreadable) const;

    std::vector<PackedBitset> collude(const std::vector<std::size_t>& coalition) const;

    // Generate codeword for user
    void getLLCodeword(std::size_t user, std::vector<PackedBitset>& out) const;

private:
    std::size_t N_, c_, L_, n_, d_, total_length_, block_length_;
    double e_;
    std::vector<FingerprintingCode> components_;
    std::vector<std::vector<std::size_t>> hiddenCode_;
    std::mt19937_64 engine_;
    std::vector<std::vector<std::size_t>> permutations_;

    std::vector<std::vector<std::size_t>> createHiddenCode();
};

class TardosCodes{
    public:
    TardosCodes(std::size_t c, double e, std::size_t n = 0);

    static std::pair<PackedBitset, std::vector<double>> generateCodeWord(const std::vector<double>& probabilities);
    static std::vector<double> calculateProbabilities(std::size_t cs, double err, std::size_t new_words=0);
    static std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> generateCodeBook(const std::vector<double>& probabilities, std::size_t new_words);
    static std::vector<std::size_t> trace(const std::vector<std::vector<double>>& U, PackedBitset& y, double z);

    std::pair<PackedBitset, std::vector<double>> writeCodeWord();
    std::pair<std::vector<PackedBitset>, std::vector<std::vector<double>>> writeCodeBook(std::size_t new_words);
    std::vector<std::size_t> trace(PackedBitset& y);
    PackedBitset collude(const std::vector<std::size_t>& coalition) const;

    std::size_t c() const { return c_; }
    std::size_t n() const { return n_; }
    double e() const { return e_; }
    double k() const { return k_; }
    double Z() const { return Z_; }
    std::size_t getLength() const { return l_; }
    std::vector<double> getProbabilities() const { return probabilities_;}
    std::vector<PackedBitset> getCodeBook() const { return codeBook_; }
    std::vector<std::vector<double>> getU() const { return U_; }

    private:
    std::vector<double> probabilities_;
    std::vector<PackedBitset> codeBook_;
    std::vector<std::vector<double>> U_;

    std::size_t l_, c_, n_;
    double e_, k_, Z_;
};

} // namespace fingerprinting

#endif // FINGERPRINTING_CODES_HPP
