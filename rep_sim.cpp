#include <iostream>
#include <cstdlib>
#include <algorithm>

class RepProb {
public:
    RepProb(unsigned d, unsigned n) : d(d), n(n) {
        // distribution[i][j] <-> ith die, jth rep
        distribution = new double*[n];
        for (std::size_t i = 0; i < d; i++) {
            distribution[i] = new double[d];
        }

        // distribution : each column is a die. each row is a boundary of [0,1]
        // i.e., distribution 0.2, 0.6, 0.9 means that [0, 0.2] -> 1, [0.2, 0.6] -> 2, [0.6, 0.9] -> 3
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < d; j++) {
                distribution[i][j] = rand();
            }
        }

        for (std::size_t i = 0; i < n; i++) {
            std::sort(distribution, distribution + d);
        }

        // Change structure so that distribution[i][j] = P(ith die = jth rep)
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 1; j < d; j++) {
                distribution[i][j] -= distribution[i][j - 1];
            }
        }

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
            std::sort(threads, threads + n, [=](unsigned dieOne, unsigned dieTwo) {
                // Sort based on each die's P( = i)
                return distribution[dieOne][i] < distribution[dieTwo][i]; 
            });
        }

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
        for (std::size_t i = 0; i < n; i++) {
            delete[] threads[i];
        }
        delete[] threads;
    }

    unsigned d, n;
    double** distribution;
    unsigned* rrOrdering;
};


int main() {
    return 0;
}