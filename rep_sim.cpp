#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <cmath>
#include <functional>

// How many permutations to check before updating the display
#define LOCAL_DISPLAY 100
#define OPT_DISPLAY 1000

int popcount(unsigned x) { //couldn't be bothered to import? 
    int count = 0;
    while (x) {
        x &= (x - 1);  // clear the lowest set bit
        ++count;
    }
    return count;
}

std::size_t factorial(std::size_t n) {
    std::size_t res = 1;
    for (std::size_t i = 2; i <= n; i++) {
        res *= i;
    }

    return res;
}

std::ostream& displayDistribution(std::ostream& os, const double* const* distribution, unsigned n, unsigned d) {
    os << "Distribution: (row : die, column : color)\n";
    for (std::size_t i = 0; i < n; i++) {
        for (std::size_t j = 0; j < d; j++) {
            os << std::round(1000.0f * distribution[i][j]) / 1000.0f << '\t';
        }
        os << '\n';
    }

    return os;
}

std::ostream& displayFavorabilities(std::ostream& os, const double* const* distribution, unsigned n, unsigned d) {
    os << "Favorabilities:\n";
    for (std::size_t color = 0; color < d; color++) {
        double favor = 0;
        for (std::size_t die = 0; die < n; die++) {
            favor += distribution[die][color];
        }

        os << color << ": " << std::round(1000.0f * favor) / 1000.0f << '\n';
    }

    return os;
}

// #define DEBUG

class RepProb {
    struct Realization {
        const unsigned* outcome = nullptr;
        double P = 1; // P(x = ordering)
    };

    friend std::ostream& operator<<(std::ostream& os, const RepProb& repProb) {
        displayDistribution(os, repProb.distribution, repProb.n, repProb.d);
        displayFavorabilities(os, repProb.distribution, repProb.n, repProb.d);

        os << "RR:\n";
        for (std::size_t i = 0; i < repProb.n; i++) {
            os << repProb.RR[i] << ' ';
        }
        os << "\nE[RR]: " << repProb.ERR << '\n';

        if (repProb.optVerified) {
            os << "OPT:\n";
            for (std::size_t i = 0; i < repProb.n; i++) {
                os << repProb.OPT << ' ';
            }
            os << "\nE[OPT]: " << repProb.EOPT;
        } else {
            // os << "Best local optimum found:\n";
            // for (std::size_t i = 0; i < repProb.n; i++) {
            //     os << repProb.localOPT[i] << ' ';
            // }
            // os << "With expectation: " << repProb.ElocalOPT;
        }

        return os;
    }
public:
    RepProb(unsigned d, unsigned n) : d(d), n(n),
    distribution(nullptr), RR(nullptr), ERR(0), OPT(new unsigned[n]),
    EOPT(n), optVerified(false), localOPT(new unsigned[n]), ElocalOPT(n),
    realizations(nullptr), numRealizations(0) {
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
        // numRealizations = pow(d, n);
        // realizations = genRealizations();

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
        ERR = expected_markov(RR);
    }

    ~RepProb() {
        delete[] RR;
        for (std::size_t i = 0; i < n; i++) {
            delete[] distribution[i];
        }
        delete[] distribution;
        delete[] OPT;
        delete[] localOPT;
        if (realizations) {
            delete[] realizations;
        }
    }

    void initDistribution() {
        numRealizations = pow(d, n);
        realizations = genRealizations();
    }

    // Returns P( |pred = true), as well as P(pred = true)
    std::tuple<double**, double> conditionalDistribution(const std::function<bool(const unsigned*)>& pred) const {
        double pPred = 0; // Probability of the predicate being true
        double** condDist = new double*[n];
        for (std::size_t i = 0; i < n; i++) {
            condDist[i] = new double[d];
            
            // Initially zero
            for (std::size_t j = 0; j < d; j++) {
                condDist[i][j] = 0.0f;
            }
        }

        // Calculate individual PMFs of the dice
        for (std::size_t i = 0; i < numRealizations; i++) {
            if (pred(realizations[i].outcome)) {
                pPred += realizations[i].P;

                for (std::size_t j = 0; j < n; j++) {
                    // jth die
                    unsigned whichRep = realizations[i].outcome[j];
                    condDist[j][whichRep] += realizations[i].P;
                }
            }
        }

        // Final scaling
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j < d; j++) {
                condDist[i][j] /= pPred;
            }
        }

        return {condDist, pPred};
    }

    void RRvsOPT() {
        exhaustiveSearch();
    }

    void RRvsLocalOPT() {
        // Not sure why C++ makes RNG so weird
        std::random_device rd;
        std::mt19937 gen(rd());

        std::size_t numTries = 1; // How many random starting points to try
        unsigned* start = new unsigned[n];
        for (std::size_t i = 0; i < n; i++) {
            start[i] = i;
        }

        for (std::size_t i = 0; i < numTries; i++) {
            std::shuffle(start, start + n, gen); // Random starting point

            double thisCost = localSearch(start); // Modifies start in-place

            // If we found a better local minimum, replace localOPT
            if (thisCost < ElocalOPT) {
                ElocalOPT = thisCost;
                for (std::size_t i = 0; i < n; i++) {
                    localOPT[i] = start[i];
                }
            }
        }

        // Cleanup
        delete[] start;
    }

    double getERR() const {
        return ERR;
    }

    double getElocalOPT() const {
        return ElocalOPT;
    }

    void testExpected() const {
        unsigned* ordering = new unsigned[n];
        for (unsigned i = 0; i < n; ++i) {
            ordering[i] = i;
        }
        std::cout << "expected(): " << expected(ordering) << std::endl; 
        std::cout << "expected() markov: " << expected_markov(ordering) << std::endl;

        delete[] ordering;
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

    double expected_markov(const unsigned* ordering) const {
        double E = 0;
        std::size_t states = 1 << d;  // 2^d using bitshift

        double* markov = new double[states]();
        bool* alive = new bool[states]();
        // Initialize
        for (std::size_t i = 0; i < states; ++i) {
            alive[i] = true;
            markov[i] = 0.0; 
        }
        alive[states - 1] = false;  // last state 1 * [d] is initially dead
        markov[0] = 1.0;  // start with full weight at state 0
        unsigned cost = 0;
    
        for (std::size_t i = 0; i < n; ++i) {
            cost ++; 
            unsigned dice = ordering[i];
    
            for (int j = states -1; j >= 0; --j) {  // reverse iteration: j from states-1 to 0

                if (!alive[j]) continue; //if not alive we have 0 anyway 
    
                double newself = 0.0;
    
                for (std::size_t k = 0; k < d; ++k) {
                    unsigned destination = j | (1u << k);  // set bit k
                    double chance = markov[j] * distribution[dice][k];
    
                    if (destination != j) {
                        if (alive[destination]) {
                            markov[destination] += chance;
                        } else {
                            E += chance * cost;
                        }
                    } else {
                        newself += chance;
                    }
                }
                markov[j] = newself;
            }
    
            // Begin killing nodes
            std::size_t left = n - i - 1; //iteration left to travel 

            if (left < d - 1) {
                for (int j = states -1; j >= 0; j--) {
                    //on the last iteration or the state is has not enough 1 
                    if (i == n - 1 || popcount(j) < left) {
                        E += markov[j] * cost;
                        alive[j] = false;
                        markov[j] = 0;
                    }
                }
            }
        }
    
        delete[] markov;
        delete[] alive;
    
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
        // Must start sorted for std::next_permutation to be exhaustive
        unsigned* ordering = new unsigned[n];
        for (std::size_t i = 0; i < n; i++) {
            ordering[i] = i;
        }

        std::size_t numChecked = 0;
        std::size_t totalPermuations = factorial(n);

        do {
            double thisCost = expected_markov(ordering);
            if (thisCost < EOPT) {
                // Copy into OPT
                for (std::size_t i = 0; i < n; i++) {
                    OPT[i] = ordering[i];
                }

                // Note cost
                EOPT = thisCost;
            }

            numChecked++;
            if (0 == numChecked % OPT_DISPLAY) {
                std::cout << "Checked " << numChecked 
                << " permutations out of " << totalPermuations << " total. "
                << std::round(1000.0f * double(numChecked) / double(totalPermuations)) / 10.0f
                << "% done.\r" << std::flush;
            }
        } while (std::next_permutation(ordering, ordering + n));
        optVerified = true;

        // Cleanup
        delete[] ordering;

        // Display what we've found
        std::cout << "OPT:\n";
        for (std::size_t i = 0; i < n; i++) {
            std::cout << OPT[i] << ' ';
        }
        std::cout << "\nE[OPT]:\n";
        std::cout << EOPT << '\n';
    }

    // Searches in the neighborhood of start for a local minimum
    double localSearch(unsigned* start) {
        // Track the best expectation we've found so far
        double expToBeat = expected_markov(start);

        // If we don't find any better neighbours, the loop will terminate
        std::size_t numChecked = 0;

        // Swapping start[bestSwapFrom] and start[bestSwapTo]
        // will be the best neighbour of start
        std::size_t bestSwapFrom = n, bestSwapTo = n;
        do {
            bestSwapFrom = bestSwapTo = n; // Reset
            for (std::size_t i = 0; i < n; i++) {
                for (std::size_t j = 0; j < n; j++ ) {
                    if (j == i) { continue; }

                    // Try swapping i with j
                    unsigned aux = start[i];
                    start[i] = start[j];
                    start[j] = aux;

                    // Check if this is better
                    double thisExp = expected_markov(start);
                    if (thisExp < expToBeat) {
                        expToBeat = thisExp;
                        bestSwapFrom = i;
                        bestSwapTo = j;
                    }

                    // Undo the change
                    start[j] = start[i];
                    start[i] = aux;
                }
            }

            // Make the change we found to be best
            if (bestSwapFrom != n) {
                unsigned aux = start[bestSwapFrom];
                start[bestSwapFrom] = start[bestSwapTo];
                start[bestSwapTo] = aux;
            }

            numChecked++;
            if (0 == numChecked % LOCAL_DISPLAY) {
                std::cout << "Checked " << numChecked << " permutations. Best: " << expToBeat << '\r';
            }
        } while (bestSwapFrom != n);

        // Let the caller know how well we did
        return expToBeat;
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

    unsigned* OPT;
    double EOPT;
    bool optVerified;

    unsigned* localOPT;
    double ElocalOPT;

    Realization* realizations;
    std::size_t numRealizations;
};


int main() {
    unsigned d = 3;

    for (unsigned n = 12; n <= 12; n++) {
        RepProb repProb(d, n);
        
        // Evaluates the representative problem
        auto f = [d, n](const unsigned* outcome){
            unsigned needCount = d;

            bool* seen = new bool[d];
            for (std::size_t i = 0; i < d; i++) {
                seen[i] = false;
            }

            unsigned testsLeft = n;
            for (std::size_t i = 0; i < n; i++) {
                if (! seen[outcome[i]]) {
                    needCount--;
                    seen[outcome[i]] = true;
                }

                testsLeft--;
                if (needCount > testsLeft) {
                    delete[] seen;
                    return false;
                }
            }

            delete[] seen;
            return 0 == needCount;
        };

        repProb.initDistribution();
        auto [dist, prob] = repProb.conditionalDistribution(f);
    
        std::cout << repProb << std::endl;
        
        std::cout << "CONDITIONED ON f = 1:\n";
        std::cout << "P(f = 1): " << std::round(prob * 1000.0f) / 1000.0f << '\n';
        displayDistribution(std::cout, dist, n, d);
        displayFavorabilities(std::cout, dist, n, d);
        std::cout << std::endl;


        //repProb.testexpected();

        // repProb.RRvsLocalOPT();
        // std::cout << repProb << '\n';
        // std::cout << "Approximation Factor: " << repProb.getERR() / repProb.getElocalOPT() << '\n';
        // std::cout << "d = " << d << " ; n = " << n << '\n';
        // std::cout << "==========================" << std::endl;
    }

    return 0;
}