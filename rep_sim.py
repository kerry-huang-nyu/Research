import numpy as np
import itertools as it
from tqdm import tqdm

class RepProb:
    def __init__(self, d, n):
        self.d = d
        self.n = n

    def check_opt_vs_rr(self):
        # distribution : each column is a die. each row is a boundary of [0,1]
        # i.e., distribution 0.2, 0.6, 0.9 means that [0, 0.2] -> 1, [0.2, 0.6] -> 2, [0.6, 0.9] -> 3
        self.distribution = np.random.rand(self.d, self.n)
        self.distribution.sort(axis=0)

        opt, Eopt = self.opt_search()
        rr_ordering = self.rr()
        ERR = self.expected(rr_ordering)

        factor = ERR / Eopt
        if factor > self.d:
            print("\n")
            print('dist: ', self.distribution)
            print('opt: ', opt)
            print('E[opt]: ', Eopt)
            print('RR: ', rr_ordering)
            print('E[RR]: ', ERR)
            print("\n")

    def generate_outcomes_for_stopping(self, order):
        # replacing it.product(range(self.d), repeat=len(order))
        # after this, we no longer need to track the set() in expected since we are doing all the work in generate orders 
        outcomes = set()

        def backtrack(path, seen):
            if len(seen) == self.d:
                outcomes.add(tuple(path))
                return

            for val in range(self.d):
                if len(path) >= len(order):  # Prevent overflow
                    return
                backtrack(path + [val], seen | {val}) #unions seen and {val} through a new set 

        backtrack([], set())
        return outcomes

    def expected(self, order):  
        # Abandon all hope ye who enter here
        # each "bit" could be {0, 1, ... d-1} 
        # why go for efficiency when you can go for a d^n algo? 
        # stop early once we find all representatives by the order 

        total_expected = 0.0

        #combinations = set() #all combinations of dice seen. 
        #if [1, 0, 2, 2, 2] has been seen [1, 0, 2] will be stored to prevent [1, 0, 2, 1, 0] from being calculated 
        # All combinations of outcomes (one per test)

        for outcome in self.generate_outcomes_for_stopping(order): #generate vector len n with values = {1 - d}
            steps = 0
            prob = 1.0 #compute probability until we see a 1 or 0 certificate 

            for i in range(len(outcome)):
                #dice is the nth dice we are pulling {0 ... n - 1}
                #val is the value on that dice we generated  {0 ... d-1}
                #sampling the ith value because it is easier to slide outcome that way 
                dice = order[i]
                val = outcome[i]

                steps += 1
                chance = self.distribution[0, dice] if val == 0 else self.distribution[val, dice] - self.distribution[val - 1, dice]
                total = self.distribution[self.d-1, dice]
                prob *= chance / total

            total_expected += prob * steps

        return total_expected

    def opt_search(self):
        current_min = float('inf')
        current_best = None

        # Brute force through all orderings
        for order in tqdm(it.permutations(range(self.n)), total=self.d ** self.n, desc="searching_for_opt"): #wait is this supposed to be a length of 5? 
            current_exp = self.expected(order)
            if current_exp < current_min:
                current_min = current_exp
                current_best = order
        
        return current_best, current_min

    # Round-robin ordering optimized for c in [d]
    def rr_for_c(self, c):
        assert 1 <= c <= self.d
        # Copy to preserve original
        prob_of_c = self.distribution[c - 1, :].copy()

        # Necessary because of the structure of distribution
        if c > 1:
            prob_of_c -= self.distribution[c - 2, :]
        
        # Descending order
        return prob_of_c.argsort()[::-1]

    def rr(self):
        threads = []
        for c in range(1, self.d + 1):
            threads.append(self.rr_for_c(c))
        
        turn = 0
        tested = set()
        rr_ordering = []
        while turn < self.n:
            for thread in threads:
                # Don't add the same test twice
                if thread[turn] not in tested:
                    rr_ordering.append(thread[turn])
                    tested.add(thread[turn])
            turn += 1
        
        return tuple(rr_ordering)


if __name__ == '__main__':
    d = 3
    n = 10
    rep_prob = RepProb(d, n)
    for _ in tqdm(range(1_00), desc="Running simulations"):
        rep_prob.check_opt_vs_rr()