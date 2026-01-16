// Copyright (c) 2026 David Brand
// SPDX-License-Identifier: MIT

#include "TTT.hpp"
#include "fingerprinting_codes.hpp"
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

void testFingerprinting_codes() {
    using namespace fingerprinting;

    std::vector<std::size_t> Ns = {1000, 10000};
    std::vector<double> Es = {0.1, 0.001};
    std::vector<float> Deltas = {0.0f, 0.1f, 0.5f};
    std::vector<std::size_t> coalition_sizes = {50, 100};

    std::mt19937_64 rng(std::random_device{}());

    for (auto wordcount : Ns) {
        for (auto p : Es) {
            for(auto delta : Deltas){
            int i = 0;
            while (i < coalition_sizes.size()) {
                auto t = coalition_sizes[i];
                i++;
                if(wordcount == 1000 && (t == 50 || delta == 0.0f || delta == 0.1f)){
                    continue;
                }
                std::vector<std::size_t> colluders;
                colluders.reserve(t);

                std::uniform_int_distribution<std::size_t> dist(0, wordcount - 1);
                for (std::size_t g = 0; g < t; ++g) {
                    colluders.push_back(dist(rng));
                }

                std::cout << "N: " << wordcount
                          << " c: " << t
                          << " e: " << p 
                          << " delta: " << delta
                          << '\n';

                try{
                auto creation_start = std::chrono::high_resolution_clock::now();
                LogLengthCodes a(wordcount, t, p, delta);
                //TardosCodes a(t, p);
                //a.writeCodeBook(wordcount);
                auto creation_stop = std::chrono::high_resolution_clock::now();
                std::cout << "code length: " << a.getLength() << '\n';
                
                auto wordCreation_start = std::chrono::high_resolution_clock::now();
                std::vector<PackedBitset> word;
                a.getLLCodeword(0, word);
                //a.generateCodeWord(a.getProbabilities());
                auto wordCreation_stop = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> wordCreation_time = wordCreation_stop - wordCreation_start;
                std::cout << "single word creation time: " << wordCreation_time.count() << '\n';

                std::chrono::duration<double> creation_time = creation_stop - creation_start;
                std::cout << "creation time: " << creation_time.count() << '\n';
                
                auto collusion_start = std::chrono::high_resolution_clock::now();
                //a.writeCodeBook(colluders.size());
                auto illegal_word = a.collude(colluders);
                auto collusion_stop = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> collusion_time = collusion_stop - collusion_start;
                std::cout << "collusion time: " << collusion_time.count() << '\n';

                auto tracing_start = std::chrono::high_resolution_clock::now();
                std::size_t guilty_user = a.trace(illegal_word, std::vector<PackedBitset>(illegal_word.size(), PackedBitset(a.componentLength())));
                //std::vector<std::size_t> accused = a.trace(illegal_word);
                auto tracing_stop = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> tracing_time = tracing_stop - tracing_start;
                //std::cout << "coalition: [";
                bool inC = false;
                for (std::size_t i = 0; i < colluders.size(); ++i) {
                    //std::cout << colluders[i] << (i + 1 < colluders.size() ? ", " : "");
                    if(colluders[i] == guilty_user) inC = true;
                    //for(auto acc : accused){
                        //if(colluders[i] == acc) inC = true;
                    //}
                }
                //std::cout << "]\n";
                //std::cout << "accused: [";
                //for (std::size_t i = 0; i < accused.size(); ++i) {
                    //std::cout << accused[i] << (i + 1 < accused.size() ? ", " : "");
                //}
                //std::cout << (inC ? " (in coalition)" : " (not in coalition)") << "]\n";
                std::cout << "guilty user: " << guilty_user << (inC ? " (in coalition)" : " (not in coalition)") << '\n';
                std::cout << "tracing time: " << tracing_time.count() << "\n\n";
                } catch (const std::exception& e) {
                    std::cout << "Error during processing: " << e.what() << "\n\n";
            }
        }
    }
    }
    }
}