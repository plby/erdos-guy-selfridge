import numpy as np
import matplotlib.pyplot as plt
import math


# A variant of Figure 3 of the original paper.



def sum( m ):

    sum=0
    for i in range(1,m+1):
        if i % 6 == 1 or i % 6 == 5:
            sum += 3/i
        sum -= 5/(5*i-1)
    return sum

def plot():
    base = range(1,100)
    exact = [sum(m) for m in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, exact, label='Discrepancy' )
    plt.title('$\\sum_{n \\leq K\'; (n,6)=1} \\frac{3}{n} - \\sum_{m \\leq K\'} \\frac{1}{n-0.2}$')
    plt.xlabel('$K\'$')
    plt.legend()
    plt.grid(True)
    plt.show()

plot()
