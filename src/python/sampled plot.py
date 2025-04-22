import numpy as np
import matplotlib.pyplot as plt
import math

base = []
t_lower = []
t_upper = []

with open('../../data/t_large_n_bounds_2.txt', 'r') as f:
    next(f)  # Skip the header line
    for line in f:
        parts = line.split()
        if len(parts) == 3:
            base.append(int(parts[0]))
            t_lower.append(int(parts[1]))
            t_upper.append(int(parts[2]))

        
def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def logfac(N):
    sum = 0
    for i in range(1, N+1):
        sum += math.log(i)
    return sum

# Test if  \sum_{p > \frac{t}{\sqrt{t}+1}} \left\lfloor \frac{N}{p} \right\rfloor \log \left( \frac{p}{t} \left\lceil \frac{t}{p} \right\rceil \right) > \log N! - N \log t

def criterion( N, t ):
    sum = 0
    threshold = logfac(N) - N * math.log(t) + 0.0001 
    tfail = t/(math.floor(math.sqrt(t)))
    for p in range(1, N + 1):
        if p > tfail:
            if is_prime(p):
                sum += math.floor(N/p) * math.log((p/t)*math.ceil(t/p))
                if sum > threshold:
                    return True
    return False

def best_t( N, init ):
    t = init
    while not criterion(N, t):
        t += 5
    while criterion(N,t):
        t -= 1
    print(f"Testing {N}: needed to increment by {t-init}") 
    return t


def plot():
    lb = [t_lower[i]/base[i] for i in range(len(base))]
    ub = [t_upper[i]/base[i] for i in range(len(base))]
    lin = [best_t(base[i], t_lower[i])/base[i] for i in range(len(base))]
    third = [1/3 for _ in base]
    c0 = 0.3044190
    asym = [1/math.e - c0 / math.log(N) for N in base]
    c1 = 0.75554808
    asym2 = [1/math.e - c0 / math.log(N) - c1 / math.log(N)**2 for N in base]
    plt.figure(figsize=[8, 6])
    plt.plot(base, lb, label='Linear programming (lower)', marker='.',linestyle='None', markersize=3 )
    plt.plot(base, ub, label='Linear programming (upper)', marker='.', linestyle='None', markersize=3) 
    plt.plot(base, lin, label='Lemma 5.1', marker='.', linestyle='None', markersize=3) 
    plt.plot(base, third, label='$1/3$', linestyle='--')
    plt.plot(base, asym, label='$1/e - c_0/\\log N$', linestyle='--')
    plt.plot(base, asym2, label='$1/e - c_0/\\log N - c_1/ \\log^2 N$', linestyle='--')
    plt.title('Approximations to $t[N]/N$')
    plt.xlabel('$N$')
    plt.ylim(0.31,0.34)
    plt.legend()
    plt.grid(True)
    plt.show()

plot()