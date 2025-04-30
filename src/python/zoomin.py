import numpy as np
import matplotlib.pyplot as plt
import math
import re

# Plot the linear algebra surplus near 43632


# data using the smoothfac method

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

# data using LP and floor function
base2 = []
t_lower2 = []
t_upper2 = []
t_lower_ip = []

with open('../../data/lp_bounds_40000-45000.txt', 'r') as f:
    next(f)
    for line in f:
        nums = re.findall(r"-?\d+(?:\.\d+)?", line)
        a, b, c, d, e = int(nums[0]), int(nums[1]), int(nums[2]), int(nums[3]), float(nums[4])
        base2.append(a)
        t_lower2.append(c)
        t_lower_ip.append(d)
        t_upper2.append(e)
        if d < a and a >= 43632:
            print(f"Ad hoc IP lower bound {d} for {a} is below the floor {c} at N={a}!")
        


def plot():
    lb = [t_lower[i]-base[i] for i in range(len(base))]
    ub = [t_upper[i]-base[i] for i in range(len(base))]
    lb2 = [t_lower2[i]-base2[i] for i in range(len(base2))]
    ub2 = [t_upper2[i]-base2[i] for i in range(len(base2))]
    lb_ip = [t_lower_ip[i]-base2[i] for i in range(len(base2))]
    plt.figure(figsize=[8, 6])
#    plt.plot(base, ub, label='Smoothfac (upper)', marker='.', markersize=3, color='green') 
    plt.plot(base2, ub2, label='Linear programming (upper)', marker='.', markersize=3, color='red') 
    plt.plot(base, lb, label='Smooth factorization (lower)', marker='.',linestyle='None', markersize=3, color='purple' )
    plt.plot(base2, lb_ip, label='Impatient IP (lower)', marker='.',linestyle='None', markersize=3, color='green' )
    plt.plot(base2, lb2, label='LP floor (lower)', marker='.',linestyle='None', markersize=3, color='blue' )
#    plt.plot(base, [0 for _ in base], label='$0$',  color='purple')
    plt.title('Surplus factors $M(N)-N$ in $N/3$-admissible factorization of N!')
    plt.xlabel('$N$')
    plt.ylabel('$M(N)-N$')
    plt.ylim(-12,3)
    plt.legend()
    plt.grid(True)
    plt.show()

plot()