#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <cmath>

#define DEBUG

class RepProb {
    friend std::ostream& operator<<(std::ostream& os, const RepProb& repProb) {
        os << "Distribution: (row : die, column : color)\n";
        for (std::size_t i = 0; i < repProb.n; i++) {
            for (std::size_t j = 0; j < repProb.d; j++) {
                os << std::round(1000.0f * repProb.distribution[i][j]) / 1000.0f << '\t';
            }
            os << '\n';
        }

        os << "RR Ordering:\n";
        for (std::size_t i = 0; i < repProb.n; i++) {
            os << repProb.rrOrdering[i] << ' ';
        }

        return os;
    }
public:
    RepProb(unsigned d, unsigned n) : d(d), n(n) {
        // distribution[i][j] <-> ith die, jth rep
        distribution = new double*[n];
        for (std::size_t i = 0; i < n; i++) {
            distribution[i] = new double[d];
        }

        // https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0f, 1.0f);

        // distribution : each column is a die. each row is a boundary of [0,1]
        // i.e., distribution 0.2, 0.6, 0.9 means that [0, 0.2] -> 1, [0.2, 0.6] -> 2, [0.6, 0.9] -> 3
        for (std::size_t i = 0; i < n; i++) {
            // Fill except last number, since we only need d - 1 partitions
            for (std::size_t j = 0; j < d - 1; j++) {
                distribution[i][j] = dis(gen);
            }
        }

        for (std::size_t i = 0; i < n; i++) {
            // Sort, except last number
            // 0 through (d - 1) are the partitions
            std::sort(distribution[i], distribution[i] + d - 1);
        }

        // Change structure so that distribution[i][j] = P(ith die = jth rep)
        for (std::size_t i = 0; i < n; i++) {
            distribution[i][d - 1] = 1.0f - distribution[i][d - 2];

            // Iterate backwards so we don't mess up future iterations
            // Skip last element, since it's already correct
            for (std::size_t j = d - 2; j > 0; j--) {
                distribution[i][j] -= distribution[i][j - 1];
            }
        }

        // Re-sort
        for (std::size_t i = 0; i < n; i++) {
            std::sort(distribution[i], distribution[i] + d);
        }

        // Initialize RR ordering
        genRR();
    }

    ~RepProb() {
        delete[] rrOrdering;
        for (std::size_t i = 0; i < n; i++) {
            delete[] distribution[i];
        }
        delete[] distribution;
    }

private:
    // Initializes rrOrdering
    void genRR() {
        // threads[i][j] = jth die in ordering optimized for ith rep
        unsigned** threads = new unsigned*[d];
        for (std::size_t i = 0; i < d; i++) {
            threads[i] = new unsigned[n];

            // Starting in the initial order
            for (std::size_t j = 0; j < n; j++) {
                threads[i][j] = j;
            }

            // Swap indices until it's the ordering the thread wants
            std::sort(threads[i], threads[i] + n, [this, i](unsigned dieOne, unsigned dieTwo) {
                // Sort based on each die's P( = i)
                // Reversed
                return distribution[dieOne][i] > distribution[dieTwo][i]; 
            });
        }

        #ifdef DEBUG
        for (std::size_t i = 0; i < d; i++) {
            std::cout << "S" << i << '\n';
            for (std::size_t j = 0; j < n; j++) {
                std::cout << threads[i][j] << ' ';
            }
            std::cout << '\n';
        }
        #endif

        rrOrdering = new unsigned[n];       // rrOrdering[i] = which die goes (i + 1)th ?
        std::size_t turn = 0;               // what position in the threads are we ?
        std::size_t rrOrderingIndex = 0;    // what position in the rrOrdering are we filling ?

        // used[i] = have we used die i ?
        bool* used = new bool[n];
        for (std::size_t i = 0; i < n; i++) {
            used[i] = false;
        }

        while (turn < n && rrOrderingIndex < n) {
            // Pull a die from each thread
            for (std::size_t whichThread = 0; whichThread < d; whichThread++) {
                unsigned candidateDie = threads[whichThread][turn];

                // Ensure the die hasn't been rolled yet
                if (! used[candidateDie]) {
                    rrOrdering[rrOrderingIndex++] = candidateDie;
                    used[candidateDie] = true;
                }
            }
            turn++;
        }

        // Cleanup
        for (std::size_t i = 0; i < d; i++) {
            delete[] threads[i];
        }
        delete[] threads;
    }

    unsigned d, n;
    double** distribution;
    unsigned* rrOrdering;
};


int main() {
    RepProb repProb(3, 10);

    std::cout << repProb << std::endl;

    return 0;
}