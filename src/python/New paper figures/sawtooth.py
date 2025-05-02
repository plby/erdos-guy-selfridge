import numpy as np
import matplotlib.pyplot as plt
import math

# Plot numbers that are coprime to 6

def coprime_6(n):
    if n % 6 == 1 or n % 6 == 5:
        return True
    return False

def sawtooth(x):
    sum = 0
    for i in range(0, math.floor(x)+1):
        if coprime_6(i):
            sum += 1
    return sum - x / 3

def plot():
    A = 30
    N = 100
    base = [i/N for i in range(0,A*N+1)]
    saw = [sawtooth(x) for x in base]
    lower = [-2/3 for _ in base]
    upper = [2/3 for _ in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, saw, label='$\\sum_{k \\leq x}^* 1 - \\frac{x}{3}$' )
    plt.plot(base, upper, label='$2/3$' )
    plt.plot(base, lower, label='$-2/3$' )
    plt.title('Normalized count of numbers coprime to 6')
    plt.xlabel('$x$')
    plt.legend()
    plt.grid(True)
    plt.show()

plot()
