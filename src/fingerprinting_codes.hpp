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
    FingerprintingCode(std::size_t n, double e, std::mt19937_64& engine, std::vector<std::size_t>& permutation, std::size_t d_override = 0);

    std::size_t n() const { return n_; }
    double e() const { return e_; }
    std::size_t d() const { return d_; }
    std::size_t l() const { return l_; }

    // Generate codeword for user
    void getCodeword(std::size_t user, PackedBitset& out) const;

    //returns indices of guilty users
    std::vector<std::size_t> trace(const PackedBitset& x) const;

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
    LogLengthCodes(std::size_t N, std::size_t c, double e);

    std::size_t N() const { return N_; }
    std::size_t c() const { return c_; }
    std::size_t L() const { return L_; }
    std::size_t n() const { return n_; }
    std::size_t d() const { return d_; }

    //return index of most likely guilty word
    std::size_t trace(const std::vector<PackedBitset>& x) const;

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

} // namespace fingerprinting

#endif // FINGERPRINTING_CODES_HPP
