import numpy as np
import matplotlib.pyplot as plt
import math



A = 110
K = 280
N = 10 ** 13
t = N / 3
L = 4.5

def kappa(L):
    if L >= 341.34:
        return math.log(9/8)
    if L >= 40.5:
        return math.log(32/27)
    if L >= 4.5:
        return math.log(4/3)
    return math.log(3/2)

kappa_K = kappa(K)
kappa_L = kappa(L)
kappa_5 = kappa(5)

print(f"kappa({K}) = {kappa_K}")
print(f"kappa({L}) = {kappa_L}")
print(f"kappa(5) = {kappa_5}")

sigma = 3*N/(t*A)
delta = math.log(3)-1

# An effective error term in the prime number theorem

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

def nu(p,m):
    count = 0
    n = m
    while n % p == 0:
        count += 1
        n //= p
    return count

def coprime_6(n):
    return n%6 == 1 or n%6 == 5

# An effective upper bound for pi(N)

def pi_upper(N):
    assert N > 1, "Error: this bound is only valid for N > 1"
    return N/math.log(N) + 1.2762 * N / (math.log(N) ** 2)


def pixy_upper(y,x):
    assert N > 1423, "Error: this bound is only valid for N > 1423"
    return (x-y)/math.log(y) + 2 * E(x) / math.log(y)

def pixy_lower(y,x):
    assert N > 1423, "Error: this bound is only valid for N > 1423"
    return (1-2/math.sqrt(y))*(x-y)/math.log(y) - 2 * E(x) / math.log(y)

# An upper bound for Z_p

def Z(p):
    sum = 0
    for m in range(K+1, math.floor(K*(1+sigma))+1):
        if coprime_6(m):
            sum += nu(p,m) * pixy_upper(t/K, t*(1+sigma)/m)
    return A*sum


# The integral of \lfloor N/x \rfloor from t/m to t(m-1), which can be computed explicitly when t=N/3

def nt_integ(m):
    sum = 0
    for i in range(3):
        if 3*m - i - 1 > 0:
            sum += (3*m-i-1) * (N/(3*m-i-1) - N/(3*m-i))
    return sum

def npsum_upper(p1):
    sum = 0
    for m in range(1,K+1):
        sm = N
        if m > 1:
            sm = t/(m-1)
        a = nt_integ(m) / math.log(t/m) + 2*t*E(sm)/(m*N/math.log(t/m))
        sum += nu(p1,m) * a
    return sum

def npsum_lower(p1):
    sum = 0
    for m in range(1,K+1):
        sm = N
        if m > 1:
            sm = t/(m-1)
        a = (1-2/math.sqrt(sm)) * nt_integ(m) / math.log(sm) - 2*t*E(sm)/(m*N/math.log(t/m))
        sum += nu(p1,m) * a
    return sum

def yp_upper(p1):
    sum = 0
    for m in range(1,K+1):
        if coprime_6(m):
            sum += A * nu(p1,m) * pixy_upper(t/m, t*(1+sigma)/m)
    return sum - npsum_lower(p1)

def yp_lower(p1):
    sum = 0
    for m in range(1,K+1):
        if coprime_6(m):
            sum += A * nu(p1,m) * pixy_lower(t/m, t*(1+sigma)/m)
    return sum - npsum_upper(p1)

print(f"Raw budget: {N * delta}")

budget = N * delta - 3 * math.log(N) / 2 - kappa_L * math.log(12) * N / (2 * math.log(t))

print(f"Budget: {budget} ( {budget / (N * delta) * 100:.2f}% of raw budget)")

def excess():
    excess_1 = 9*N/(2*A) + 4
    excess_2 = (N / math.log(N/(3*K))) * 0.9201
    excess_3 = (N / math.log(N/(3*K))) * 2044 * E(N) / N
    return excess_1 + excess_2 + excess_3

excess = excess()
print(f"Excess: {excess} ({excess / budget * 100:.2f}% of budget)")


print(f"log(N/3K) = {math.log(N/(3*K))}; log(N) = {math.log(N)}")

X1_1 = ((4*A+3)/3) * pi_upper(t/K)
X1_2 = ((4*A+3)/3) * pi_upper(math.sqrt(N)) * math.log(N) / math.log(5)
X1_3 = A * (pi(K*(1+sigma)) - pi(K)) * pixy_upper(t/K, (t*(1+sigma))/K)

X1_4 = 0
for p in range(4,K+1):
    if is_prime(p):
        X1_4 += Z(p) * math.log(p) / math.log(math.sqrt(t/K))

X1_5 = 0
for p in range(4,K+1):
    if is_prime(p):
#        norm = N / ((p-1)*math.log(N))
#        print(f"Y_{p} in [{yp_lower(p)/norm}, {yp_upper(p)/norm}]; npsum({p}) in [{npsum_lower(p)/norm}, {npsum_upper(p)/norm}]")
        X1_5 += max(yp_upper(p),0) * math.log(p) / math.log(math.sqrt(t/K))
        
X1 = X1_1 + X1_2 + X1_3 + X1_4 + X1_5

print(f"First component of X1: {X1_1/N}; {X1_1 / X1 * 100:.2f}% of X1")
print(f"Second component of X1: {X1_2/N}; {X1_2 / X1 * 100:.2f}% of X1")
print(f"Third component of X1: {X1_3/N}; {X1_3 / X1 * 100:.2f}% of X1")
print(f"Fourth component of X1: {X1_4/N}; {X1_4 / X1 * 100:.2f}% of X1")
print(f"Fifth component of X1: {X1_5/N}; {X1_5 / X1 * 100:.2f}% of X1")


X2_1 = (4*A+3)/3 * (pi(K)-2) * (math.log(N) / math.log(5))

X2_2 = 0
for p in range(4,K+1):
    if is_prime(p):
        X2_2 += max(-yp_lower(p),0)

X2 = X2_1 + X2_2

print(f"Excess: {excess/N}")

print(f"X1: {X1/N}")

print(f"X2: {X2/N}")

print(f"budget: {budget/N}")

print(f"Excess / budget: {excess / budget * 100:.2f}%")

print(f"X1 / budget: {X1 * kappa_K / budget * 100:.2f}%")

print(f"X2 / budget: {X2 * kappa_5 / budget * 100:.2f}%")

print(f"Remaining: {(budget - excess - X1 * kappa_K - X2 * kappa_5) / budget * 100:.2f}%")

def volume():
    volume1 = N * math.log(3) / 2 - math.log(N) - ((math.log(2*L)+kappa_L) / (math.log(t) - math.log(2*L))) * N * math.log(2)
    volume2 = N * math.log(2)  - math.log(N) - ((math.log(3*L)+kappa_L) / (math.log(t) - math.log(3*L))) * N * math.log(3)/2
    return min(volume1, volume2)
    
volume = volume()

X1_volume = (math.log(math.sqrt(t*K))+kappa_K)*(X1+2)
X2_volume = (math.log(5)+kappa_5)*X2

print(f"Volume: {volume / N:.2f}")
print(f"X1 / volume: { X1_volume / volume * 100:.2f}%")

print(f"X2 / volume: { X2_volume / volume * 100:.2f}%")

print(f"Remaining: {(volume - X1_volume - X2_volume) / volume * 100:.2f}%")

