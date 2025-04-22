import numpy as np
import matplotlib.pyplot as plt
import math

base = range(1, 10001)
values=[0 for _ in base]

def read_file():
    # Open the file for reading
    with open("..\\..\\data\\a034258-lpsolve.txt", "r") as file:
        # Read the file line by line
        for line in file:
            # Split the line into parts (assuming numbers are separated by whitespace)
            parts = line.split()
            # Ensure there are at least two numbers on the line
            if len(parts) >= 2:
                values[int(parts[0])-1] = int(parts[1])



def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

primes = [is_prime(n) for n in base]

def logfac(N):
    sum = 0
    for i in range(1, N+1):
        sum += math.log(i)
    return sum

# Test if  \sum_{p > \frac{t}{\sqrt{t}+1}} \left\lfloor \frac{N}{p} \right\rfloor \log \left( \frac{p}{t} \left\lceil \frac{t}{p} \right\rceil \right) > \log N! - N \log t

def criterion( N, t ):
    sum = 0
    for p in range(1, N + 1):
        if primes[p-1] and p > t/(math.floor(math.sqrt(t))):
            sum += math.floor(N/p) * math.log((p/t)*math.ceil(t/p))
    if sum > logfac(N) - N * math.log(t) + 0.0001:
        return True


def best_t( N, init ):
    for t in range(init, N+1):
        if criterion(N,t):
            return t-1
    return N



def plot2():
    read_file()

    base2 = range(80, 10001)
    sparse_base = [10*n for n in range(8, 1001)]
    sparse_upper = [best_t(N, values[N-1])/N for N in sparse_base]

    exact = [values[N-1]/N for N in base2]
    asym = [1/math.e - 0.3044019010/math.log(N) for N in base2]
    asym2 = [1/math.e - 0.30440119010/math.log(N) - 0.75554808/math.log(N)**2 for N in base2]
    e_inv = [1/math.e for N in base2]
    third = [1/3 for N in base2]

    plt.figure(figsize=(8, 6))
    plt.plot(base2, exact, label='$t(N)/N$ (exact)' )
    plt.plot(base2, asym, label='$1/e-c_0/\\log N$' )
    plt.plot(base2, asym2, label='$1/e-c_0/\\log N$-c_1/\\log^2 N' )
    plt.plot(base2, e_inv, linestyle="--", label='$1/e$' )
    plt.plot(base2, third, linestyle="--", label='$1/3$' )
    plt.plot(sparse_base, sparse_upper, label='Lemma 5.1' )
    plt.title('$t(N)/N$')
    plt.xlabel('$N$')
    plt.ylim(0.27,0.37)
    plt.legend()
    plt.grid(True)
    plt.show()

plot2()
