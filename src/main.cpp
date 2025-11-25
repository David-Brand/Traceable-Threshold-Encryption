#include "btbf.hpp"

#include <iostream>
#include <iomanip>
#include <algorithm>

static void printKey(const btbf::BTBF::SymKey& k, const char* label)
{
    std::cout << label << ": 0x";
    for (auto b : k) {
        std::cout << std::hex << std::setw(2) << std::setfill('0')
                  << static_cast<int>(b);
    }
    std::cout << std::dec << "\n";
}

int main()
{
    try {
        btbf::BTBF::init();

        // Example parameters
        int n   = 5;
        int t   = 3;
        int ell = 2; // two positions, j=1,2

        auto keygenOut = btbf::BTBF::keygen(n, t, ell);

        int j = 1; // position to encapsulate to

        // Encapsulate
        auto [k_enc, ct] = btbf::BTBF::encaps(keygenOut.pk, j);

        // Let parties 1,2,4 contribute decryption shares using their LEFT keys at position j
        std::vector<int> J = {1, 2, 4}; // indices of parties (1-based)

        std::vector<btbf::BTBF::GT> shares;
        shares.reserve(J.size());
        for (int idx : J) {
            const auto& party = keygenOut.parties[idx - 1];
            const auto& sk_left_j = party.left[j - 1]; // left key (b=0) for position j
            auto di = btbf::BTBF::decShare(sk_left_j, ct.c0);
            shares.push_back(di);
        }

        auto k_dec = btbf::BTBF::combine(J, shares);

        printKey(k_enc, "Encapsulated key");
        printKey(k_dec, "Recovered   key");

        bool equal = std::equal(k_enc.begin(), k_enc.end(), k_dec.begin());
        std::cout << "Keys match? " << (equal ? "YES" : "NO") << "\n";

        return equal ? 0 : 1;
    }
    catch (const std::exception& ex) {
        std::cerr << "Exception: " << ex.what() << "\n";
        return 1;
    }
}
