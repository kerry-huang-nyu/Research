import numpy as np
import itertools as it

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
            print('dist: ', self.distribution)
            print('opt: ', opt)
            print('E[opt]: ', Eopt)
            print('RR: ', rr_ordering)
            print('E[RR]: ', ERR)

    def expected(self, order):
        # Abandon all hope ye who enter here
        pass

    def opt_search(self):
        current_min = float('inf')
        current_best = None

        # Brute force through all orderings
        for order in it.permutations(range(len(self.distribution))):
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
                if not tested[thread[turn]]:
                    rr_ordering.append(thread[turn])
                    tested.add(thread[turn])
            turn += 1
        
        return tuple(rr_ordering)


if __name__ == '__main__':
    d = 3
    n = 5
    rep_prob = RepProb(d, n)
    for _ in range(1_000):
        rep_prob.check_opt_vs_rr()