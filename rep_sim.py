import numpy as np
import itertools as it
from tqdm import tqdm
import math

class RepProb:
    def __init__(self, d, n):
        self.d = d
        self.n = n
        self.outcomes = self.generate_outcomes_for_stopping()

    def check_opt_vs_rr(self):
        # distribution : each column is a die. each row is a boundary of [0,1]
        # i.e., distribution 0.2, 0.6, 0.9 means that [0, 0.2] -> 1, [0.2, 0.6] -> 2, [0.6, 0.9] -> 3
        self.distribution = np.random.rand(self.d, self.n)
        self.distribution.sort(axis=0)

        rr_ordering = self.rr()
        ERR = self.expected(rr_ordering)
        print('RR: ', [ int(i) for i in rr_ordering ])
        print('E[RR]: ', ERR)
        opt, Eopt = self.opt_search()

        factor = ERR / Eopt
        if factor > self.d:
            print("\n")
            print('dist: ', self.distribution)
            print('opt: ', opt)
            print('E[opt]: ', Eopt)
            print('RR: ', rr_ordering)
            print('E[RR]: ', ERR)
            print("\n")

    def generate_outcomes_for_stopping(self):
        # replacing it.product(range(self.d), repeat=len(order)) generating d**n
        # this will only generate sequences within lenth n
        outcomes = set()

        def backtrack(path, seen):
            if len(seen) == self.d or len(path) >= self.n: #prevent overflow
                outcomes.add(tuple(path))
                return
            
            need = self.d - len(seen)
            rolls_left = self.n - len(path)

            if need <= rolls_left: #I will only keep checking if it is still possible 
                for val in range(self.d): 
                    backtrack(path + [val], seen | {val}) #unions seen and {val} through a new set 
            else: #terminate the sequence 
                outcomes.add(tuple(path))

        backtrack([], set())
        return outcomes

    def expected(self, order):  
        # Abandon all hope ye who enter here

        total_expected = 0.0

        for outcome in self.outcomes: #generate vector len n with values = {1 - d}
            # if len(outcome) > len(order): #skip the current iteration the order is less than the outcome vector 
            #     continue

            steps = len(outcome)
            prob = 1.0 #compute probability until we see a 1 or 0 certificate 

            #dice is the nth dice we are pulling {0 ... n - 1}
            #val is the value on that dice we generated  {0 ... d-1}
            #sampling the ith value because it is easier to slide outcome that way 

            for i in range(len(outcome)):
                dice = order[i]
                val = outcome[i]
                chance = self.distribution[0, dice] if val == 0 else self.distribution[val, dice] - self.distribution[val - 1, dice]
                total = self.distribution[self.d-1, dice]
                prob *= chance / total

            total_expected += prob * steps

        return total_expected

    def opt_search(self):
        current_min = float('inf')
        current_best = None

        total_perms = math.factorial(self.n)
        # Brute force through all orderings
        for i, order in enumerate(it.permutations(range(self.n))): #wait is this supposed to be a length of 5? 
            if i % 1000 == 0:
                print(f'Current min E[S]: {current_min}')
                print(f'Achieved by ordering {current_best}')
                print(f"[status] Checked {i}/{total_perms} permutations ({(i / total_perms * 100):.2f}%)", end='\r')
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
    n = 5
    rep_prob = RepProb(d, n)
    for _ in tqdm(range(100), desc="Running simulations"):
        rep_prob.check_opt_vs_rr()