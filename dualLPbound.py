import numpy as np
import itertools as it

def getQ(n, d):
    return 0

def u(S, a, d, Q):
    n = len(a)

def uS(S, a, i, d, Q):
    Si = S.copy(); Si[i] = True
    return u(Si, a, d, Q) - u(S, a, d, Q)

def expression(S, a, d, Q):
    numerator = sum(
        [ uS(S, a, i, d) for i, selected in enumerate(S) if not selected ]
    )
    return numerator / (Q - u(S, a, d, Q))

def all_subsets(a, d, Q):
    n = len(a)
    S = np.zeros(n, dtype=bool)
    max_so_far = -1
    while True:
        max_so_far = max(max_so_far, expression(S, a, d, Q))
        i = 0
        while i < n and S[i]:
            S[i] = False; i += 1
        if n == i:
            return max_so_far
        S[i] = True

def bound(n, d):
    a = np.array([1] * n, dtype=int)
    max_so_far = -1
    Q = getQ(n, d)
    while 1:
        max_so_far = max(max_so_far, all_subsets(a, d, Q))
        
        i = 1
        while i < n + 1 and d == a[-i]:
            a[-i] = 1; i += 1
        if n + 1 == i:
            return max_so_far
        a[-i] += 1

        print(a)

bound(3, 4)
all_subsets([1, 2, 3])