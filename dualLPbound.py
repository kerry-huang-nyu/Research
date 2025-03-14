import numpy as np
import itertools as it

class RepProb:
    __slots__ = (
        'd', 'n', 'Q', 'subQs',
        'S', 'a', 'maxS', 'maxa',
        'current_u', 'max_so_far'
    )

    def __init__(self, d, n):
        assert d <= n # problem should be feasible
        self.d = d
        self.n = n
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

    def u(self):
        num_tests = len(list(filter(bool, self.S)))
        unique = len({ e for e in it.compress(self.a, self.S) })
        subtractor = self.d - unique # 1-certificate
        for i in range(2, self.d - 1):
            unique_val = min(unique, i) # capped at goal value
            tests_val = min(self.n - self.d + 1 + i, num_tests) # capped at goal value
            subtractor *= self.subQs[i - 1] - (unique_val + tests_val)
        return self.Q - subtractor

    def uS(self, i):
        assert not self.S[i]

        self.S[i] = True
        uSi = self.u()
        self.S[i] = False

        assert uSi - self.current_u >= 0
        return uSi - self.current_u

    def expression(self):
        numerator = sum(
            [ self.uS(i) for i in range(len(self.S)) if not self.S[i] ]
        )
        assert self.Q - self.current_u >= 0
        return numerator / (self.Q - self.current_u)

    def all_subsets(self):
        self.S = np.zeros(self.n, dtype=bool)
        while True:
            self.current_u = self.u()
            if self.current_u != self.Q:
                new_value = self.expression()
                if new_value > self.max_so_far:
                    self.max_so_far = new_value
                    # save current S and a
                    self.maxS = self.S.copy()
                    self.maxa = self.a.copy()

            i = 0
            while i < self.n and self.S[i]:
                self.S[i] = False; i += 1
            if self.n == i:
                return
            self.S[i] = True

    def bound(self):
        while True:
            self.all_subsets()
            
            i = 1
            while i < self.n + 1 and self.d == self.a[-i]:
                self.a[-i] = 1; i += 1
            if self.n + 1 == i:
                return self.max_so_far, self.maxS, self.maxa
            self.a[-i] += 1

if __name__ == '__main__':
    # seems to be bounded at n - d + 1 ?
    for d, n in it.product(range(3, 6), range(7, 10)):
        print(f'd={d}, n={n}')
        bound, S, a = RepProb(d, n).bound()
        print(bound, S, a, sep='\n')
        print('---------------------')