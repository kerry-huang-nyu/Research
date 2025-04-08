import numpy as np
import itertools as it

class RepProb:
    __slots__ = (
        'd', 'n', 'Q', 'subQs',
        'S', 'a', 'maxS', 'maxa',
        'current_u', 'max_so_far'
    )

    def __init__(self, d, n):
        assert d <= n # problem must be feasible
        self.d, self.n = d, n
        self.initQs()
        self.a = np.ones(self.n, dtype=int)
        self.current_u = 0
        self.max_so_far = float('-inf')

    def initQs(self):
        # [0] is the 1-certificate
        # [1, d - 1] are the 0-certificates with goal values [2, d - 2]
        self.subQs = np.empty(self.d - 2)
        self.Q = self.subQs[0] = self.d # 1-certificate
        # TODO: this could probably be done better
        for i in range(1, self.d - 2): # [2, d - 2]
            # seen (i + 1) unique reps AND performed (n - d + 1) + (i + 1) tests
            thisQ = (i + 1) + (self.n - self.d + 1 + (i + 1))
            self.subQs[i] = thisQ; self.Q *= thisQ

    # returns u(S, a) for current S and a
    def u(self):
        num_tests = len(list(filter(bool, self.S)))
        num_reps = len({ e for e in it.compress(self.a, self.S) })
        # subtractor will end as the product of all the OR terms
        subtractor = self.d - num_reps # 1-certificate
        for i in range(2, self.d - 1):
            # both below are capped at goal value
            unique_val = min(num_reps, i)
            tests_val = min(self.n - self.d + 1 + i, num_tests)
            subtractor *= self.subQs[i - 1] - (unique_val + tests_val)
        return self.Q - subtractor

    # returns the marginal value of x_(i + 1)
    def uS(self, i):
        assert not self.S[i]

        self.S[i] = True
        uSi = self.u()
        self.S[i] = False

        assert uSi - self.current_u >= 0
        return uSi - self.current_u

    # calculates the fraction in the dual LP bound
    def expression(self):
        numerator = sum(
            [ self.uS(i) for i in range(self.n) if not self.S[i] ]
        )

        # this method shouldn't be called if u(S, a) = Q
        assert self.Q - self.current_u > 0
        return numerator / (self.Q - self.current_u)

    # finds the maximum over all subsets of a given assignment a
    def all_subsets(self):
        self.S = np.zeros(self.n, dtype=bool)
        while True:
            # update and cache current u(S, a)
            self.current_u = self.u()

            # only consider if Q - u(S, a) != 0
            if self.current_u != self.Q:
                new_value = self.expression()
                if new_value > self.max_so_far:
                    self.max_so_far = new_value
                    # save current S and a
                    self.maxS = self.S.copy()
                    self.maxa = self.a.copy()

            # increment to next subset
            i = 0
            while i < self.n and self.S[i]:
                self.S[i] = False; i += 1
            if self.n == i:
                return
            self.S[i] = True

    # finds the dual LP bound for the specified d and n
    def bound(self):
        while True:
            self.all_subsets()
            
            # increment to next assignment
            i = 1
            while i < self.n + 1 and self.d == self.a[-i]:
                self.a[-i] = 1; i += 1
            if self.n + 1 == i:
                return self.max_so_far, self.maxS, self.maxa
            self.a[-i] += 1

if __name__ == '__main__':
    # seems to be bounded at n - d + 1 ?
    n_domain = range(7, 11)
    d_domain = range(3, 8)
    print('d\tn', end='')
    for i in n_domain:
        print(f'\t{i}', end='')
    print()
    for d in d_domain:
        print(f'{d}\t', end='')
        for n in n_domain:
            bound, S, a = RepProb(d, n).bound()
            print(f'\t{bound}', end='')
        print()