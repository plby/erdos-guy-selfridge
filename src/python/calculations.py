import numpy as np
import matplotlib.pyplot as plt
import math


# The code here verifies the calculations in Section 7.  The bounds are monotone in N, so if the criterion is verified for one value of N, it automatically holds for larger N.



# An effective error term in the prime number theorem, see (C.9)
def E(N):
    return 0.95 * math.sqrt(N) + 2.39 * 10**(-8) * N

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

def pi_upper(N):
    assert N > 1, "Error: this bound is only valid for N > 1"
    return N/math.log(N) + 1.2762 * N / (math.log(N) ** 2)

# Upper bound for pi(x)-pi(y), see (C.5)
def pixy_upper(y,x):
    assert N > 1423, "Error: this bound is only valid for N > 1423"
    return (x-y)/math.log(y) + 2 * E(x) / math.log(y)

# Lower bound for pi(x)-pi(y), see (C.6)
def pixy_lower(y,x):
    assert N > 1423, "Error: this bound is only valid for N > 1423"
    return (1-2/math.sqrt(y))*(x-y)/math.log(y) - 2 * E(x) / math.log(y)

# See Table 2
def kappa(L):
    if L >= 341.34:
        return math.log(9/8)
    if L >= 40.5:
        return math.log(32/27)
    if L >= 4.5:
        return math.log(4/3)
    if L >= 4/3:
        return math.log(3/2)
    return math.log(2)

def fancy_kappa(L):
    return (2/math.log(3) - 2/math.log(12))*math.log(12*L) + (2/math.log(3))*kappa(L)

def delta(t,N):
    return math.log(N/t) - 1

# see (6.3)
def sigma(t,N,A):
    return 3*N/(t*A)

# Some minor terms appearing in (6.23)
def minor_delta_terms(L,t,N):
    return kappa(L) * math.log(12) / (2 * math.log(t)) + 3 * math.log(N) / (2*N)

# (6.16)
def alpha_eval(N,t,L):
    alpha1 = (math.log(3*L) + kappa(L))/(math.log(t) - math.log(3*L)) * math.log(3)/(2*math.log(2)) + math.log(N) / (N*math.log(2)) + 1/N
    alpha2 = (math.log(2*L) + kappa(L))/(math.log(t) - math.log(2*L)) * 2 * math.log(2) / math.log(3) + 2 * math.log(N) / (N * math.log(3)) + 2/N
    print(f"alpha1={alpha1}, alpha2={alpha2}")
    return max(alpha1, alpha2)

# equation (6.36)    
def excess_bound(t,N,A,K):
    sum = 3*N/(2*t*A)
    sum += 4/N
    sum += 0.9201 / math.log(t/K)
    sum += 2044 * E(N) / (N * math.log(t/K))
    return sum

# equation (6.32)
def Z_bound(t,N,A,K,p):
    sum = 0
    for m in range(K+1, math.floor(K*(1+sigma(t,N,A)))+1):
        if coprime_6(m):
            sum += nu(p,m) * pixy_upper(t/K, t*(1+sigma(t,N,A))/m) / N
    return sum

# The first term in (6.37)
def Y1p_first(t,N,A,K):
    sum = 0
    for p in range(4, K+1):
        if is_prime(p):
            sum += (math.log(N)/math.log(p) + 1) * math.log(p) 
    sum *= (4*A + 3) / (3 * N * math.log(t/K**2))
    print(f"First term in bound on Y1+: {sum}")
    return sum

# The first term in (6.38)
def Y1m_first(t,N,A,K):
    sum = 0
    for p in range(4, K+1):
        if is_prime(p):
            sum += (math.log(N)/math.log(p) + 1)
    sum *= (4*A + 3) / (3 * N)
    print(f"First term in bound on Y1-: {sum}")
    return sum

# The first term in (6.39)
def Y2pm_first(t,N,A,K):
    sum = 0
    sum += pi_upper(t/K)
    sum += pi_upper(math.sqrt(N)) * (math.log(N)/math.log(K) + 1)
    sum *= (4*A + 3) / (3 * N)
    print(f"First term in bound on Y2+-: {sum}")
    return sum

# The expression in (6.41)
# assumes t = N/3
def tinyprimes_bound(t,N,K):
    sum2 = 0
    sum3 = 0
    for m in range(1, K+1):
        x = (3*m-1) * pixy_upper(t/(3*m), t/(3*m-1))
        x += (3*m-2) * pixy_upper(t/(3*m-1), t/(3*m-2))
        if 3*m-3 > 0:
            x += (3*m-3) * pixy_upper(t/(3*m-2), t/(3*m-3))
        sum2 += x * nu(2,m) / N
        sum3 += x * 2*nu(3,m) / N
    return max(sum2,sum3)

# upper bound for (6.41)
# assumes t = N/3
def Wp_upper(t,N,A,K,p):
    sum = 0
    for m in range(1, K+1):
        if coprime_6(m):
            sum += (A/N) * nu(p,m) * pixy_upper(t/m, t*(1+sigma(t,N,A))/m)
    
    for m in range(1, K+1):
        x = (3*m-1) * pixy_lower(t/(3*m), t/(3*m-1))
        x += (3*m-2) * pixy_lower(t/(3*m-1), t/(3*m-2))
        if 3*m-3 > 0:
            x += (3*m-3) * pixy_lower(t/(3*m-2), t/(3*m-3))
        sum -= nu(p,m) * x / N
    return sum

# lower bound for (6.41)
# assumes t = N/3
def Wp_lower(t,N,A,K,p):
    sum = 0
    for m in range(1, K+1):
        if coprime_6(m):
            sum += (A/N) * nu(p,m) * pixy_lower(t/m, t*(1+sigma(t,N,A))/m)
    
    for m in range(1, K+1):
        x = (3*m-1) * pixy_upper(t/(3*m), t/(3*m-1))
        x += (3*m-2) * pixy_upper(t/(3*m-1), t/(3*m-2))
        if 3*m-3 > 0:
            x += (3*m-3) * pixy_upper(t/(3*m-2), t/(3*m-3))
        sum -= nu(p,m) * x / N
    if sum < 0:
        print(f"Possible negative value of W_{p}")
    return sum

# (6.37)        
def Y1p_bound(t,N,A,K):
    sum = Y1p_first(t,N,A,K)
    zsum=0
    wsum=0
    for p in range(4, K+1):
        if is_prime(p):
            zsum += Z_bound(t,N,A,K,p) * math.log(p) / math.log(t/K**2)
            wsum += max(0, Wp_upper(t,N,A,K,p)) * math.log(p) / math.log(t/K**2)
    print(f"Contribution of Z_p terms to Y1+: {zsum}")
    print(f"Contribution of W_p terms to Y1+: {wsum}")
    return sum+zsum+wsum

# (6.38)
def Y1m_bound(t,N,A,K):
    sum = Y1m_first(t,N,A,K)
    for p in range(4, K+1):
        if is_prime(p):
            sum += max(0, -Wp_lower(t,N,A,K,p))
    return sum

# (6.39)
def Y2pm_bound(t,N,A,K):
    sum = Y2pm_first(t,N,A,K)
    for p in range(K+1, math.floor(K*(1+sigma(t,N,A)))+1):
        if is_prime(p):
            sum += (A/N) * pixy_upper(t/K, t*(1+sigma(t,N,A))/p)
    return sum

# Check if (6.24), (6.25) hold.
def evaluate(N, A, K, L):
    t = N/3
    print("(N,A,K,L,t,sigma) = ", N, A, K, L, t, sigma(t,N,A))
    assert t > 9*L, "Error: t must be greater than 9*L"
    assert K**2*(1+sigma(t,N,A)) < t, "Error: K^2*(1+sigma) must be less than t"
    assert K*math.sqrt(N) < t, "Error: K*sqrt(N) must be less than t"
    assert sigma(t,N,A) < 1, "Error: sigma must be less than 1"
    assert L >= 4.5, "Error: L must be at least 4.5"
    assert K >= L, "Error: K must be at least L"


    Y1p = Y1p_bound(t, N, A, K)
    Y1m = Y1m_bound(t, N, A, K)
    Y2pm = Y2pm_bound(t, N, A, K)

    d = delta(t, N)

# This term is like 1/A; to make it smaller, increase A
    delta1 = excess_bound(t, N, A, K)

# This term is like log K / log N; to make it smaller, decrease K or increase N
    delta2 = kappa(4.5) * Y1p

# This term is negligible
    delta3 = kappa(4.5) * Y1m

# This term is like A / K log N; to make it smaller, decrease A, increase K, or increase N
    delta4 = kappa(4.5) * Y2pm

# this term is like kappa(L) / log N; to make it smaller, increase N or increase L
    delta5 = minor_delta_terms(L, t, N)

    total = delta1 + delta2 + delta3 + delta4 + delta5
    slack = d - total

    print(f"delta budget: {d}")
    print(f"Excess term: {delta1} ({delta1 / d * 100:.2f}% of delta)")
    print(f"Y1p term: {delta2} ({delta2 / d * 100:.2f}% of delta)")
    print(f"Y1m term: {delta3} ({delta3 / d * 100:.2f}% of delta)")
    print(f"Y2pm term: {delta4} ({delta4 / d * 100:.2f}% of delta)")
    print(f"Minor delta terms: {delta5} ({delta5 / d * 100:.2f}% of delta)")
    print(f"Total left-hand side of delta equation: {total} ({total / d * 100:.2f}% of delta)")

# Q decreases if L increases
    Q = 1 - alpha_eval(N, t, L)

# This term is like log K (log log K) / log N; to make it smaller, decrease K or increase N
    Q1 = tinyprimes_bound(t, N, K)

# This term is like log K (log log K) / log N; to make it smaller, decrease K or increase N
    Q2 = (2/math.log(12)) * (math.log(K**2) + fancy_kappa(4.5)) * Y1p

# This term is negligible
    Q3 = (2/math.log(12)) * (math.log(K) + fancy_kappa(4.5)) * Y1m

# This term is like A/K; to make it smaller, decrease A or increase K
    Q4 = (2/math.log(12)) * (math.log(t/K) + fancy_kappa(K)) * Y2pm

# A negligible term
    Q5 = (2/math.log(12)) * (math.log(t) + fancy_kappa(K)) / N

    total2 = Q1 + Q2 + Q3 + Q4 + Q5
    slack2 = Q - total2

    print(f"Tiny-prime budget: {Q}")
    print(f"Direct term: {Q1} ({Q1 / Q * 100:.2f}% of 1-\\alpha)")
    print(f"Y1p term: {Q2} ({Q2 / Q * 100:.2f}% of 1-alpha)")
    print(f"Y1m term: {Q3} ({Q3 / Q * 100:.2f}% of 1-alpha)")
    print(f"Y2pm term: {Q4} ({Q4 / Q * 100:.2f}% of 1-alpha)")
    print(f"Negligible term: {Q5} ({Q5 / Q * 100:.2f}% of 1-alpha)")
    print(f"Total left-hand side of tiny prime equation: {total2} ({total2 / Q * 100:.2f}% of 1-\\alpha)")
    assert slack > 0, "Over budget for delta"
    assert slack2 > 0, "Over budget for tiny primes"
    print(f"Guy-Selfridge conjecture t(N) >= N/3 successfully verified for N >= {N}")


A = 175
K = 340
N = 10 ** 12
L = 4.5

evaluate(N, A, K, L)
