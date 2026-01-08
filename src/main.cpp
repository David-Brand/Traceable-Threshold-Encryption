#include "btbf.hpp"
#include "TTT.hpp"

#include <chrono>

void testTTT_BTBF(int n, int t, int lambda, double e){
    using namespace ttt;
    TTT_BTBF E;

    std::cout << "TTT_BTBF with n: " << n << " t: " << t << " lambda: " << lambda << " e: " << e << "\n";

    auto keygen_start = std::chrono::high_resolution_clock::now();
    auto keygenOut = E.keygen(n, t, lambda, e);
    auto keygen_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> keygen_duration = keygen_end - keygen_start;

    std::vector<int> J(t);
    std::iota(J.begin(), J.end(), 0);

    auto enc_start = std::chrono::high_resolution_clock::now();
    auto cip = E.enc(keygenOut.public_key);
    auto enc_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> enc_duration = enc_end - enc_start;

    auto dec_start = std::chrono::high_resolution_clock::now();
    std::vector<GT> shares;
    shares.reserve(J.size());
    for (int idx : J) {
        shares.push_back(E.dec(keygenOut.parties[idx], idx, cip.second));
    }
    auto k_rec = E.combine(J, shares);
    auto dec_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dec_duration = dec_end - dec_start;

    auto fingerprint_length = E.fingerPrintingCode().getLength();

    std::cout << "Fingerprint Length:  " << fingerprint_length << "\n";
    std::cout << "Encryption Time:     " << enc_duration.count() << "\n";
    std::cout << "Decryption Time:     " << dec_duration.count() << "\n";

    auto D = [&J, &keygenOut, &E](const Ciphertext& c, const SymKey& k) {
        // Simulated decoder box oracle
        std::vector<GT> shares;
        shares.reserve(J.size());
        for (int idx : J) {
            shares.push_back(E.dec(keygenOut.parties[idx], idx, c));
        }
        auto k_rec = E.combine(J, shares);
        return k == k_rec;
    };

    auto trace_start = std::chrono::high_resolution_clock::now();
    std::vector<std::size_t> accused = E.trace(keygenOut.public_key, keygenOut.tracing_key, D);
    auto trace_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> trace_duration = trace_end - trace_start;

    auto total_time = keygen_duration + trace_duration;

    std::cout << "Key Generation Time: " << keygen_duration.count() << "\n";
    std::cout << "Tracing Time:        " << trace_duration.count() << "\n";
    std::cout << "Total Time:          " << total_time.count() << "\n";
}

void testTTT_BTDDH(int n, int t, int lambda, double e){
    using namespace ttt;
    TTT_BTDDH E;

    std::cout << "TTT_BTDDH with n: " << n << " t: " << t << " lambda: " << lambda << " e: " << e << "\n";
    
    auto keygen_start = std::chrono::high_resolution_clock::now();
    auto keygenOut = E.keygen(n, t, lambda, e);
    auto keygen_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> keygen_duration = keygen_end - keygen_start;

    std::vector<int> J(t);
    std::iota(J.begin(), J.end(), 0);

    auto enc_start = std::chrono::high_resolution_clock::now();
    auto cip = E.enc(keygenOut.public_key);
    auto enc_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> enc_duration = enc_end - enc_start;

    auto dec_start = std::chrono::high_resolution_clock::now();
    std::vector<G> shares;
    shares.reserve(J.size());
    for (int idx : J) {
        shares.push_back(E.dec(keygenOut.parties[idx], idx, cip.second));
    }
    auto k_rec = E.combine(J, shares);
    auto dec_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dec_duration = dec_end - dec_start;

    auto fingerprint_length = E.fingerPrintingCode().getLength();

    std::cout << "Fingerprint Length:  " << fingerprint_length << "\n";
    std::cout << "Encryption Time:     " << enc_duration.count() << "\n";
    std::cout << "Decryption Time:     " << dec_duration.count() << "\n";

    auto D = [&J, &keygenOut, &E](const Ciphertext_ddh& c, const G& k) {
        // Simulated decoder box oracle
        std::vector<G> shares;
        shares.reserve(J.size());
        for (int idx : J) {
            shares.push_back(E.dec(keygenOut.parties[idx], idx, c));
        }
        auto k_rec = E.combine(J, shares);
        return k == k_rec;
    };

    auto trace_start = std::chrono::high_resolution_clock::now();
    std::vector<std::size_t> accused = E.trace(keygenOut.public_key, keygenOut.tracing_key, D);
    auto trace_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> trace_duration = trace_end - trace_start;

    auto total_time = keygen_duration + trace_duration;

    std::cout << "Fingerprint Length:  " << fingerprint_length << "\n";
    std::cout << "Key Generation Time: " << keygen_duration.count() << "\n";
    std::cout << "Tracing Time:        " << trace_duration.count() << "\n";
    std::cout << "Total Time:          " << total_time.count() << "\n";
}

int main(){
    using namespace ttt;

    bool testBTDDH = false;
    bool testBTBF = true;

    std::mt19937_64 rng{ std::random_device{}() };
    int ns[] = {100, 200};
    int ts[] = {25, 100};

    int lambdas[] = {4};
    double es[] = {0.99, 0.5};

    for(int n : ns){
        for(int t : ts){
            for(int lambda : lambdas){
                for(double e : es){
                    if(testBTDDH){
                        testTTT_BTDDH(n, t, lambda, e);
                    }
                    if(testBTBF){
                        testTTT_BTBF(n, t, lambda, e);
                    }
                }
            }
        }
    }
    return 0;
}
