#ifndef FINGERPRINTING_CODES_HPP
#define FINGERPRINTING_CODES_HPP

#include <vector>
#include <cstddef>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>

namespace fingerprinting {

using Bitset = std::vector<bool>;

class FingerprintingCode {
public:
    FingerprintingCode(std::size_t n, double e, std::size_t d_override = 0);

    // Accessors
    std::size_t n() const { return n_; }
    double e() const { return e_; }
    std::size_t d() const { return d_; }
    std::size_t l() const { return l_; }

    void getCodeword(std::size_t user, Bitset& out) const;
    void getCodewordToSlice(std::size_t user, Bitset& out, std::size_t offset) const;

    // Operations
    std::vector<std::size_t> trace(const Bitset& x) const;
    Bitset collude(std::size_t i, std::size_t j) const;


private:
    std::size_t n_;                 // number of users
    double      e_;                 // error probability
    std::size_t d_;                 // block length
    std::size_t l_;                 // total length = d * (n-1)
    std::vector<std::size_t> permutation_;

    std::size_t weight(const Bitset& x,
                       std::size_t start,
                       std::size_t end) const; // end exclusive

    // global RNG helper
    static std::mt19937_64& rng();
};


// Log-length codes wrapper
class LogLengthCodes {
public:
    LogLengthCodes(std::size_t N, std::size_t c, double e);

    // Accessors
    std::size_t N() const { return N_; }
    std::size_t c() const { return c_; }
    std::size_t L() const { return L_; }
    std::size_t n() const { return n_; }
    std::size_t d() const { return d_; }

    // Trace pirate copy and return index of most likely guilty word
    std::size_t trace(const std::vector<Bitset> x) const;

    // Coalition collusion: given indices of colluders, produce illegal word
    std::vector<Bitset> collude(const std::vector<std::size_t>& coalition) const;

    std::vector<Bitset> getLLCodeword(std::size_t user) const;

private:
    std::size_t N_; // number of words
    std::size_t c_; // coalition size
    std::size_t L_; // number of component codes
    std::size_t n_; // n = 2c
    double      e_;
    std::size_t d_; // block length in each FingerprintingCode
    std::size_t total_length_;
    std::size_t block_length_;

    std::vector<FingerprintingCode> components_;
    std::vector<std::vector<std::size_t>> hiddenCode_;

    std::vector<std::vector<std::size_t>> createHiddenCode();
};

} // namespace fingerprinting

#endif // FINGERPRINTING_CODES
