#include "btbf.hpp"
#include <iostream>
#include <vector>

using namespace btbf;

int main() {
    BTBF::init();

    // Demo params
    const int n = 5;
    const int t = 3;
    const int ell = 2; // positions

    PublicKey pk;
    std::vector<PartySecretKey> sks;
    BTBF::KeyGen(n, t, ell, pk, sks);

    // Encapsulate to position j = 1
    Ciphertext ct;
    SharedKey kEnc;
    BTBF::Enc(pk, 1, ct, kEnc);

    // Get t decryption shares
    std::vector<int> J = {1, 2, 3};              // parties indices (1-based)
    std::vector<GT> shares;
    for (int idx : J) {
        // here only left key (b=0)
        shares.push_back(BTBF::DecShare(sks[idx-1], 1, /*b=*/0, ct));
    }

    SharedKey kDec = BTBF::Combine(J, shares);

    std::cout << "Encapsulated key: " << kEnc.hex << "\n";
    std::cout << "Decapsulated key: " << kDec.hex << "\n";
    std::cout << "Match: " << (kEnc.hex == kDec.hex ? "YES" : "NO") << std::endl;
    return 0;
}
