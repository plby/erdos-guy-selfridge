import numpy as np
import matplotlib.pyplot as plt
import math
import re

# Plot the linear algebra surplus near 43632

base = []
t_lower = []
t_upper = []

with open('../../data/crossover_bounds.txt', 'r') as f:
    for line in f:
        nums = re.findall(r"-?\d+(?:\.\d+)?", line)
        a, b, c = int(nums[0]), int(nums[1]), float(nums[2])
        base.append(a)
        t_lower.append(b)
        t_upper.append(c)
        
for i in range(len(base)):
    if base[i] % 3 == 1:
      if base[i] > 43632:
        if t_lower[i] < base[i]:
            print(f"At N={base[i]} (N/3={base[i]/3}), a defect of {t_lower[i]-base[i]}")  


def plot():
    lb = [t_lower[i]-base[i] for i in range(len(base))]
    ub = [t_upper[i]-base[i] for i in range(len(base))]
    plt.figure(figsize=[8, 6])
    plt.plot(base, lb, label='Linear programming (lower)', marker='.',linestyle='None', markersize=3, color='blue' )
    plt.plot(base, ub, label='Linear programming (upper)', marker='.', markersize=3, color='green') 
    plt.plot(base, [0 for _ in base], label='$0$',  color='purple')
    plt.title('Surplus factors in $N/3$-admissible factorization of N!')
    plt.xlabel('$N$')
#    plt.ylim(0.31,0.34)
    plt.legend()
    plt.grid(True)
    plt.show()

plot()