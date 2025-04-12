import numpy as np
import matplotlib.pyplot as plt
import math


# A variant of Figure 3 of the original paper.



def sum( m ):

    sum=0
    for i in range(1,m+1):
        if i % 6 == 1 or i % 6 == 5:
            sum += 3/i
        sum -= 5 * math.log((5*i)/(5*i-1))
    return sum

def plot():
    base = range(1,100)
    exact = [sum(m) for m in base]

    plt.figure(figsize=(8, 6))
    plt.plot(base, exact, label='Discrepancy' )
    plt.title('$\\sum_{m \\leq M; (m,6)=1} 3/m - \\sum_{m \\leq M} 5 \\log \\frac{5m}{5m-1} $')
    plt.xlabel('$M$')
    plt.legend()
    plt.grid(True)
    plt.show()

plot()
