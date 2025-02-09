import numpy as np

# Setup
n = 1000
coins = np.random.rand(n)
coins = coins.sort()

# @param    bias is the probability of heads
# @return   True = Heads, False = Tails
def flip_biased(bias: float) -> bool:
    return np.random.rand() < bias # < vs <= is arbitrary here

first = np.random.randint(0, n)
flips = 1
if flip_biased(coins[first]):
    # Got heads, looking for tails
    for i, coin in enumerate(coins):
        if i == first: continue # Flipped this coin already
        flips += 1
        if not flip_biased(coin): break # If we get tails, we're done
else:
    # Got tails, looking for heads
    for i, coin in enumerate(coins[::-1]):
        if i == first: continue
        flips += 1
        if flip_biased(coin): break # If we get heads, we're done

print("n:", n)
print("coins:", coins)
print("FLIPS:", flips)