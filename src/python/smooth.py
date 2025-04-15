import numpy as np
import matplotlib.pyplot as plt
import math

# Analysis of an algorithm in which

# 1. one starts with the standard factorization $1 \times \dots \times N$
# 2. One removes all primes up to $W$ from the factors (this frees up about $N/(p-1)$ powers of p for each $p \leq W$ as the "budget").
# 3. One enumerates the $W$-smooth numbers in increasing order as $1=d_1<d_2<\dots$
# 4. One inspects the factors currently in $[N/d_k, N/d_{k+1}]$ (they occur with multiplicity $k$) and finds the first $W$-smooth $m$ that the budget permits to multiply by to move these factors to exceed $N/3$


W = 5  # has to be at least 3, may well be prime

# Limit of how much $d_k$ can reach
K = 2000

# returns true of n is W-smooth (all prime factors at most W)
def is_smooth(n):
    for d in range(2,W+1):
        while n % d == 0:
            n //= d
        if n == 1:
            return True
    return False

# the next W-smooth number greater than n
def next_smooth(n):
    while True:
        n += 1
        if is_smooth(n):
            return n
        

def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

# number of times p divides n
def nu(p,n):
    count = 0
    while n % p == 0:
        count += 1
        n //= p
    return count

primes = []
for p in range(1,W+1):
    if is_prime(p):
        primes.append(p)

density = 1
for p in primes:
    density *= (1 - 1/p)


def afford(multiplier, multiplicity, budget, k, next_k):
    for p in primes:
        if budget[p] < density * multiplicity * nu(p,multiplier) * (1/k - 1/next_k):
            return False
    return True


def calc1():

    budget = {}
    value = 0

    # At any given time, N `budget`[p] is the asymptotic estimate of how many powers of p are remaining.  N `value`` is the estimate of the logarithm of all the primes remaining, and is a measure of the total value of the budget

    for p in primes:
        budget[p] = 1/(p-1)
        value += budget[p] * math.log(p)

    # Of the N factors, `remaining` is the proportion of the factors that are still to be allocated
    remaining = 1

    k = 1
    next_k = 1
    multiplier = 1
    multiplicity = 1

    while k < K:
        next_k = next_smooth(k)
        multiplier = next_smooth(k//3)

    # greedily select the next multiplier one can afford
        while not afford(multiplier, multiplicity, budget, k, next_k):
            multiplier = next_smooth(multiplier)
            assert multiplier <= K, "Error: multiplier exceeds K"

 #   print(f"In [N/{next_k}, N/{k}], multiplier = {multiplier} and multiplicity = {multiplicity}")

        for p in primes:
            budget[p] -= density * multiplicity * nu(p,multiplier) * (1/k - 1/next_k)
        remaining -= density * multiplicity * (1/k - 1/next_k)
        value -= density * multiplicity * math.log(multiplier) * (1/k - 1/next_k)

        assert value > remaining*math.log(multiplier), "Error: budget value is not sufficient"
        k = next_k
        multiplicity += 1
        print(k, budget, remaining,value - remaining*math.log(multiplier))

A = 10

def calc2():

    budget = {}
    value = 0

    # At any given time, N `budget`[p] is the asymptotic estimate of how many powers of p are remaining.  N `value`` is the estimate of the logarithm of all the primes remaining, and is a measure of the total value of the budget

    for p in primes:
        budget[p] = 1/(p-1)
        value += budget[p] * math.log(p)
    num_smooth = 0
    remaining = 1

    for k in range(2, K+1):
        num_rough = density / (k * (k-1))
        if is_smooth(k):
            num_smooth += 1
        remaining -= num_rough * num_smooth
        for d in range(1,k):
            if is_smooth(d):
                for p in primes:
                    budget[p] -= num_rough * nu(p,d)
                value -= num_rough * math.log(d)

        print(k, remaining, budget, value)    


calc2()