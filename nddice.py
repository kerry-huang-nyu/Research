import numpy as np
import math
from itertools import product

def round_print(A: np.array):
    print(A.round(3))

def OR(a: np.float64, b: np.float64) -> np.float64:
    return 1 - (1 - a) * (1 - b)

def NOT(p: np.float64) -> np.float64:
    return 1 - p

def expected_one_die(P: np.array, n):
    d = len(P)
    if n < d: return n
    prod_all = P.prod()

    # pneed[i, j] = P(done due to rolling i on (d + j)th roll)
    phaveall = np.empty((d, n - d), dtype=np.float64)
    # pendon[i] = P(rolling exactly d + i times)
    pendon = np.empty(n - d + 1, dtype=np.float64)

    phaveall[:, 0] = prod_all
    pendon[0] = prod_all * d # phaveall[:, 0].sum()
    for j in range(1, n - d): # (d + j)th roll
        for i in range(d): # ith outcome
            phaveall[i, j] = phaveall[i, j - 1] * (1 - P[i])
        pendon[j] = phaveall[:, j].sum()
    # P(ending on nth roll) = 1 - P(ending on another roll)
    pendon[n - d] = 1 - pendon[:n - d].sum()

    return sum([ (d + i) * p for i, p in enumerate(pendon) ])

# @param    P is a matrix of probabilities such that entry i,j corresponds to
#           the probability of the jth die resulting in i
# @return   The expected number of rolls needed after ordering the dice by
#           the ordering of their corresponding columns in P
def expected(P: np.array) -> np.float64:
    d, n = P.shape
    if n < d: return n
    
    # P(Only missing i just before (d + )jth roll)
    pneed = np.empty((d, n), dtype=np.float64)
    # p_all_at_d = 
    pneed[:, 0] = []

    # P(Exactly i of face j right after roll k)
    # phave = np.empty((n, d, n + 1), dtype=np.float64)
    # phave[0, :, 0] = 1 # exactly 0 of face j right after roll 0, i.e., right before roll 1
    # phave[1:, :, 0] = 0 # can't have anything before roll 1
    # for k in range(1, n):
    #     for i, j in product(range(k + 1), range(d)):
    #         pone_less = phave[i - 1, j, k - 1] * P[k - 1, j]
    #         phad_enough = phave[i, j, k - 1] * NOT(P[k - 1, j])
    #         phave[i, j, k] = OR(pone_less, phad_enough)
    #     phave[k + 1:, :, k] = 0

    # pmissing = np.empty((d, n + 1))
    # pmissing[:, 0] = 1
    
    
    # r = np.empty(n, dtype=np.float64)
    # r[:d] = 1; r[d] = NOT(phave[1, :, d].sum())
    # for i in range(d + 1, n):

P = np.random.rand(4)
for n in range(7, 80):
    print(n, expected_one_die(P, n))