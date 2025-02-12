#include <iostream>
#include <random>

#define L 5000000
#define T 1000

double simulate_gigaflips(int num_simulations = 1) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> rng(0, L - 1);

    double mean = 0;
    for (int sim = 0; sim < num_simulations; ++sim) {
        bool booleans[L] = {};
        int flips = 0;
        int unflipped_count = L;

        while (unflipped_count > T) {
            int start = rng(gen);
            int flip_end = start + T < L ? start + T : L;

            for (int i = start; i < flip_end; ++i) {
                if (!booleans[i]) {
                    booleans[i] = true;
                    --unflipped_count;
                }
            }

            ++flips;
        }
         
        mean += flips;

    }

    return mean/num_simulations;
}

int main() {
    int num_simulations = 3;

    double mean_flips = simulate_gigaflips(num_simulations);
    std::cout << mean_flips << std::endl;

    return 0;
}

