import numpy as np
import matplotlib.pyplot as plt
import math


values = []
forced = []

def largest_prime_factor(n):
    p = 2
    if n < 1:
        return 0
    while True:
        while n % p == 0:
            n //= p
        if n == 1:
            return p
        p += 1

t = 14545
N = 43632

primes = [0 for _ in range(N + 1)]

def is_forced(n):
    p = largest_prime_factor(n)
    if p < t / math.floor(math.sqrt(t)):
        return False
    if n//p == math.ceil(t/p):
        primes[p] += 1
        return True
    else:
        return False

def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def val(p):
    m = N // p
    count = m
    while m > 0:
        m //= p
        count += m
    return count

def read_file():
    # Open the file for reading
    with open("..\\..\\..\\data\\factorizations\\43632-43632-14545.txt", "r") as file:
        # Read the file line by line
        for line in file:
            parts = line.split()
            if len(parts) == 1:
                n = int(parts[0])
                if is_forced(n):
                    forced.append(n)
                else:
                    values.append(n)

read_file()

print(len(forced), len(values))
print(f"Maximum value: {max(values)}")

for p in range(1, N + 1):
    if is_prime(p):
        if p >= t / math.floor(math.sqrt(t)):
            v = val(p)
            assert primes[p] <= v, f"Overflow at prime {p}: {primes[p]} > {v}"
            if primes[p] < v:
                print(f"Prime {p} is not entirely forced: {primes[p]} < {v}")

bins = np.linspace(14545, 16000, 100)

# Compute histograms manually to control bar placement
counts1, _ = np.histogram(values, bins)
counts2, _ = np.histogram(forced, bins)

# Midpoints of bins for placing bars
bin_centers = 0.5 * (bins[:-1] + bins[1:])
width = (bins[1] - bins[0]) / 3  # width of each bar (smaller than bin width)

# Plot
plt.bar(bin_centers + width/2, counts2, width=width, label='Factors of the form $p \\lceil t/p \\rceil$ for $p \\geq t/\\lceil \\sqrt{t} \\rceil$')
plt.bar(bin_centers - width/2, counts1, width=width, label='Other factors')

plt.title('Optimal factorization of $N!$ for $N=43632$')
plt.xlabel('Value')
plt.ylabel('Count')
plt.legend()
plt.show()
