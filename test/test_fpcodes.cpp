#include "fingerprinting_codes.hpp"
#include <iostream>
#include <random>
#include <chrono>

int main() {
    using namespace fingerprinting;

    std::vector<std::size_t> Ns = {100};
    std::vector<double> Es = {0.2};
    std::vector<std::size_t> coalition_sizes = {50};

    std::mt19937_64 rng(std::random_device{}());

    for (auto wordcount : Ns) {
        for (auto p : Es) {
            for (auto t : coalition_sizes) {
                std::vector<std::size_t> colluders;
                colluders.reserve(t);

                std::uniform_int_distribution<std::size_t> dist(0, wordcount - 1);
                for (std::size_t g = 0; g < t; ++g) {
                    colluders.push_back(dist(rng));
                }

                std::cout << "N: " << wordcount
                          << " c: " << t
                          << " e: " << p << '\n';

                auto creation_start = std::chrono::high_resolution_clock::now();
                LogLengthCodes a(wordcount, t, p);
                auto creation_stop = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> creation_time = creation_stop - creation_start;
                std::cout << "creation time: " << creation_time.count() << '\n';

                auto collusion_start = std::chrono::high_resolution_clock::now();
                std::vector<Bitset> illegal_word = a.collude(colluders);
                auto collusion_stop = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> collusion_time = collusion_stop - collusion_start;
                std::cout << "collusion time: " << collusion_time.count() << '\n';

                auto tracing_start = std::chrono::high_resolution_clock::now();
                std::size_t guilty_user = a.trace(illegal_word);
                auto tracing_stop = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> tracing_time = tracing_stop - tracing_start;
                std::cout << "coalition: [";
                for (std::size_t i = 0; i < colluders.size(); ++i) {
                    std::cout << colluders[i] << (i + 1 < colluders.size() ? ", " : "");
                }
                std::cout << "]\n";
                std::cout << "guilty user: " << guilty_user << '\n';
                std::cout << "tracing time: " << tracing_time.count() << "\n\n";
            }
        }
    }

    return 0;
}
