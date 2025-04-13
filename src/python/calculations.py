import numpy as np
import matplotlib.pyplot as plt
import math



A = 150
K = 342
N = 10 ** 11
t = N / 3
L = 342

kappa_K = math.log(9/8)
kappa_L = math.log(9/8)
kappa_5 = math.log(4/3)
sigma = 9/A
delta = math.log(3)-1

def E(N):
    E = 0.95 * math.sqrt(N)
    if N > 10**19:
        E += (math.sqrt(N) / (8*math.pi)) * math.log(N) * (math.log(N) - 3)
    if N > math.exp(45):
        E = 1.11742 * 10**(-8) * N
    return E


def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def pi(N):
    """Return the number of primes less than or equal to N."""
    count = 0
    for i in range(2, math.floor(N) + 1):
        if is_prime(i):
            count += 1
    return count

def pi_upper(N):
    return N/math.log(N) + 1.2762 * N / (math.log(N) ** 2)

print(f"Raw budget: {N * delta}")

budget = N * delta - 3 * math.log(N) / 2 - kappa_L * math.log(12) * N / (2 * math.log(t))

print(f"Budget: {budget} ( {budget / (N * delta) * 100:.2f}% of raw budget)")

excess_1 = 9*N/(2*A) + 4

excess_2 = (N / math.log(N/(3*K))) * 0.9201

excess_3 = (N / math.log(N/(3*K))) * 2044 * E(N) / N

excess = excess_1 + excess_2 + excess_3

print(f"Excess 1: {excess_1} ({excess_1 / budget * 100:.2f}% of budget)")

print(f"Excess 2: {excess_2} ({excess_2 / budget * 100:.2f}% of budget)")

print(f"Excess 3: {excess_3} ({excess_3 / budget * 100:.2f}% of budget)")

print(f"Excess: {excess} ({excess / budget * 100:.2f}% of budget)")


print(f"log(N/3K) = {math.log(N/(3*K))}; log(N) = {math.log(N)}")

X1_1 = ((4*A+3.75)/3) * pi_upper(t/K)
X1_2 = ((4*A+3.75)/3) * pi_upper(math.sqrt(N)) * math.log(N) / math.log(5)

X1_3 = A * (pi(K*(1+sigma)) - pi(K)) * (t*sigma/K) / math.log(t/K)

X1 = X1_1 + X1_2 + X1_3

print(f"First component of X1: {X1_1}; {X1_1 / X1 * 100:.2f}% of X1")
print(f"Second component of X1: {X1_2}; {X1_2 / X1 * 100:.2f}% of X1")
print(f"Third component of X1: {X1_3}; {X1_3 / X1 * 100:.2f}% of X1")




X2_1 = (4*A+3.75)/3 * (pi(K)-2) * (math.log(N) / math.log(5))
      
X2 = X2_1

print(f"X1: {X1}")

print(f"X2: {X2}")

print(f"X1 / budget: {X1 * kappa_K / budget * 100:.2f}%")

print(f"X2 / budget: {X2 * kappa_5 / budget * 100:.2f}%")


volume = N * math.log(2) + N * math.log(3) / 2

print(f"X1 / volume: {(math.log(math.sqrt(t*K))+kappa_K)*(X1+2) / volume * 100:.2f}%")

print(f"X2 / volume: {(math.log(K)+kappa_5)*X2 / volume * 100:.2f}%")
