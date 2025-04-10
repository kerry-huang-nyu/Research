#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <cmath>

// #define DEBUG

class RepProb {
    struct Realization {
        const unsigned* outcome = nullptr;
        double P = 1; // P(x = ordering)
    };

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
            os << repProb.RR[i] << ' ';
        }
        os << "\nE[RR]:\n";
        os << repProb.ERR;

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

        // Initialize realizations
        numRealizations = pow(d, n);
        realizations = genRealizations();
        #ifdef DEBUG
        std::cout << numRealizations << std::endl;
        for (std::size_t i = 0; i < numRealizations; i++) {
            for (std::size_t j = 0; j < n; j++) {
                std::cout << realizations[i].outcome[j] << ' ';
            }
            // Should flush the output buffer since there's no way
            // we get through all the realizations
            std::cout << std::endl;
        }
        #endif

        // Initialize RR ordering
        RR = genRR();
        ERR = expected(RR);
    }

    ~RepProb() {
        delete[] RR;
        for (std::size_t i = 0; i < n; i++) {
            delete[] distribution[i];
        }
        delete[] distribution;
    }

private:
    unsigned cost(const unsigned* ordering, const unsigned* realization) const {
        bool* colorsSeen = new bool[d];
        for (std::size_t i = 0; i < d; i++) {
            colorsSeen[i] = false;
        }
        unsigned uniqueSeen = 0, remainingTests = n, cost = 0;

        // Each iteration of this loop represents rolling a die
        for (std::size_t i = 0; i < n; i++) {
            cost++;
            unsigned whichDie = ordering[i];
            unsigned color = realization[whichDie];

            // If we haven't seen this color before, increment
            if (! colorsSeen[color]) {
                uniqueSeen++;
                colorsSeen[color] = true;
            }
            remainingTests--;

            // f = 1
            if (uniqueSeen == d) {
                delete[] colorsSeen;
                return cost;
            
            // Hopeless; f = 0
            } else if (remainingTests < d - uniqueSeen) {
                delete[] colorsSeen;
                return cost;
            }
        }

        delete[] colorsSeen;
        return cost;
    }

    double expected(const unsigned* ordering) const {
        double E = 0;
        for (std::size_t i = 0; i < numRealizations; i++) {
            Realization* real = realizations + i;
            E += real->P * cost(ordering, real->outcome);
        }

        return E;
    }

    // Generates and returns all d^n possible realizations
    Realization* genRealizations() const {
        // There are d^n realizations
        Realization* realizations = new Realization[numRealizations];

        // We'll modify this and use it to get the next realization
        unsigned* currentOutcome = new unsigned[n];
        
        // Start with (0,...,0)
        for (std::size_t i = 0; i < n; i++) {
            currentOutcome[i] = 0;
        }

        std::size_t curRealizationIndex = 0;
        while (true) {
            realizations[curRealizationIndex].outcome = currentOutcome;
            // Calculate probability of this particular outcome
            for (std::size_t i = 0; i < n; i++) {
                // *= P(ith die producing this color)
                realizations[curRealizationIndex].P *= distribution[i][currentOutcome[i]];
            }

            // Increment to next realization
            // First, copy the current realization
            currentOutcome = new unsigned[n];
            for (std::size_t i = 0; i < n; i++) {
                currentOutcome[i] = realizations[curRealizationIndex].outcome[i];
            }

            // i is the index of the die we're currently incrementing
            std::size_t i = 0;
            // If we increment a die and hit d, then we should carry to the next die
            // (Since really, we store values in [0, d - 1])
            while (i < n && d == ++currentOutcome[i]) {
                currentOutcome[i++] = 0;
            }
            // If we ended the prior loop because i == n, we're done
            if (i == n) {
                return realizations;
            }
            curRealizationIndex++;
        }
    }

    void exhaustiveSearch() {
        // Initially, roll the dice in reverse order
        // This makes it easier to increment to the next permutation
        unsigned* ordering = new unsigned[n];
        for (std::size_t i = 0; i < n; i++) {
            ordering[i] = n - i - 1;
        }

        unsigned* minOrdering = nullptr;
        double minE = n;
        std::size_t swapFrom = 0, swapTo = 0;
        while (true) {

            // Increment to next ordering
            if (ordering[swapTo] == n - 1) {
                swapTo++;
            }
            std::swap(ordering[swapFrom], ordering[swapTo]);
            swapTo++;
        }
    }

    // Initializes RR
    unsigned* genRR() const {
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

        unsigned* ordering = new unsigned[n];       // RR[i] = which die goes (i + 1)th ?
        std::size_t turn = 0;               // what position in the threads are we ?
        std::size_t RRIndex = 0;    // what position in the RR are we filling ?

        // used[i] = have we used die i ?
        bool* used = new bool[n];
        for (std::size_t i = 0; i < n; i++) {
            used[i] = false;
        }

        while (turn < n && RRIndex < n) {
            // Pull a die from each thread
            for (std::size_t whichThread = 0; whichThread < d; whichThread++) {
                unsigned candidateDie = threads[whichThread][turn];

                // Ensure the die hasn't been rolled yet
                if (! used[candidateDie]) {
                    ordering[RRIndex++] = candidateDie;
                    used[candidateDie] = true;
                }
            }
            turn++;
        }

        // Cleanup
        delete[] used;
        for (std::size_t i = 0; i < d; i++) {
            delete[] threads[i];
        }
        delete[] threads;

        return ordering;
    }

    unsigned d, n;
    double** distribution;
    unsigned* RR;
    double ERR;
    Realization* realizations;
    std::size_t numRealizations;
};


int main() {
    RepProb repProb(3, 7);

    std::cout << repProb << std::endl;

    return 0;
}