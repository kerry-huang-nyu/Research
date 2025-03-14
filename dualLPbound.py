import numpy as np
import itertools as it

class RP:
    __slots__ = ('d', 'n', 'Q', 'subQs', 'S', 'a')

    def __init__(self, d, n):
        self.d = d
        self.n = n
        self.initQs()
        self.a = np.ones(self.n, dtype=int)

    def initQs(self):
        self.Qs = np.empty(self.d - 2) # (d - 2) - 2 + 1 + 1
        self.Q = self.Qs[0] = self.d # 1-certificate
        # TODO: this could probably be done better
        for i in range(1, self.d - 2): # [2, d - 2]
            # seen (i + 1) unique reps AND performed n - d + 1 + (i + 1) tests
            thisQ = (i + 1) + (self.n - self.d + 1 + (i + 1))
            self.Qs[i] = thisQ; self.Q *= thisQ

    def u(self):
        subtractor = 1
        num_tests = len(list(filter(bool, self.S)))
        unique = len({ e for e in it.compress(self.a, self.S) })
        for i in range(2, self.d - 1):
            subtractor

        return self.Q - subtractor

    def uS(self, i):
        assert not self.S[i]

        u = self.u()

        self.S[i] = True
        uSi = self.u()

        self.S[i] = False
        assert uSi - u > 0
        return uSi - u

    def expression(self):
        numerator = sum(
            [ self.uS(i) for i, _ in enumerate(it.compress(self.a, not self.S)) ]
        )
        return numerator / (self.Q - self.u())

    def all_subsets(self):
        self.S = np.zeros(self.n, dtype=bool)
        max_so_far = -1
        while True:
            max_so_far = max(max_so_far, self.expression())

            i = 0
            while i < self.n and self.S[i]:
                self.S[i] = False; i += 1
            if self.n == i:
                return max_so_far
            self.S[i] = True

    def bound(self):
        max_so_far = -1
        while True:
            max_so_far = max(max_so_far, self.all_subsets())
            
            i = 1
            while i < self.n + 1 and self.d == self.a[-i]:
                self.a[-i] = 1; i += 1
            if self.n + 1 == i:
                return max_so_far
            self.a[-i] += 1