import numpy as np
import matplotlib.pyplot as plt
import math
from mpmath import *

# The code here uses interval arithmetic to compute the constants c_0, c_1.

mp.dps = 50  # Set decimal places for precision

def low(x_mpi):
    return float(x_mpi.a)

def high(x_mpi):
    return float(x_mpi.b)

# Computes the sum 1/2e \sum_{k=1}^{infinity} \log^2(1+1/k) in the formula for c0 with rigorous error bars by subtracting off the # divergent part 1/2e \sum_{k=1}^{infinity} 1/k^2 = 1/2e \zeta(2) = pi^2/12e
def c0_sum(K):
    # Use mpmath for high precision calculations
    sum_value = mpi(0,0)
    sum_value += pi**2 / 6
    for k in range(1, K + 1):
        term = log(1 + 1 / k)**2 - 1 / (k**2)
        sum_value += term

    # add uncertainty for error term, noting that log^2(1+1/k)-1/k^2 lies between -1/k^3 and 0 and summig using the integral test
    sum_value -= mpi( -1/(2*K**2), 0 )

    return sum_value / (2 * e)

# The constant term 2/e^2 - log(2)/2e in the formula for c0
def c0_const():
    return 2/e**2 - log(2)/(2*e)

# Computes the integral 1/e \int_e^infty {y} log (\lceil y/e \rceil/(y/e)) dy to high precision by evaluating the integral exactly
def c0_integral(T):
    a = e
    n = 2 # current value of \lfloor y\rfloor
    m = 2 # current value of \lceil y/e \rceil 

    integral = mpi(0,0)

    while True:
        next_n = n
        next_m = m
        if n+1 < m*e:
            next_n += 1
            next_a = n+1
        else:
            next_m += 1
            next_a = m * e
        
        if next_a > T:
            next_a = T
            terminate = True
        else:
            terminate = False

        integral += (log(m)+1) * (log(next_a) + n/next_a) - log(next_a)**2/2 - n * (log(next_a)+1) / next_a
        integral -= (log(m)+1) * (log(a) + n/a) - log(a)**2/2 - n * (log(a)+1) / a

        if terminate:
            break

        a = next_a
        n = next_n
        m = next_m

    # estimate integrand by [0,e/y^3], so tail integral is in [0, e/2T^2]
    integral += mpi(0, e/(2*T**2))
    return integral/e
        
def c0(K,T):
    return c0_sum(K) + c0_const() - c0_integral(T)
    

print(c0(100000,100000))
