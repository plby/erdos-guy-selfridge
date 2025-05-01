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

def err(x_mpi):
    return high(x_mpi) - low(x_mpi)

# Computes the sum 1/2e \sum_{k=1}^{infinity} \log^2(1+1/k) in the formula for c0 with rigorous error bars by subtracting off the # divergent part 1/2e \sum_{k=1}^{infinity} 1/k^2 = 1/2e \zeta(2) = pi^2/12e
def c0_sum(K):
    # Use mpmath for high precision calculations
    sum_value = mpi(0,0)
    sum_value += pi**2 / 6
    for k in range(1, K + 1):
        term = log(1 + 1 / k)**2 - 1 / (k**2)
        sum_value += term

    # add uncertainty for error term, noting that log^2(1+1/k)-1/k^2 lies between -1/k^3 and 0 and summing using the integral test
    sum_value -= mpi( -1/(2*K**2), 0 )

    print(f"After {K} terms, c0 sum is known to be {sum_value}; error is {err(sum_value)}")
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

    integral /= e

    print(f"On integrating up to {T}, c0 integral is known to be {integral}; error is {err(integral)}")
    return integral
        
def c0(K,T):
    return c0_sum(K) + c0_const() - c0_integral(T)



# Computes the \sum_{k=1}^{infinity} (1 + log(k+1))/2e \log^2(1+1/k) - log^3(1+1/k)/3e in the formula for c1' with rigorous error bars.  Here I don't know of a good term to subtract to speed up the convergence (it involves zeta'(2)) 
def c1p_sum(K):
    # Use mpmath for high precision calculations
    sum_value = mpi(0,0)
    for k in range(1, K + 1):
        term = (1 + log(k+1))/(2*e) * log(1 + 1 / k)**2 - log(1+1/k)**3 / (3*e)
        sum_value += term
    
    # add uncertainty for error term, noting that the summand  can be bounded above by (1 + log k)/2ek^2 + 1/6ek^3 and from below by (1+log k)/2ek^2 - (1 + log k)/2ek^3, and use the integral test

    sum_value += mpi( 1/(2*e*(K+1)) + (log(K+1)+1)/(2*e*(K+1)) - 1/(4*e*(K+1)**2) - (2*log(K+1)+1)/(8*e*(K+1)**2), 1/(2*e*K) + (log(K)+1)/(2*e*K) + 1/(12*e*K**2))
    
    print(f"After {K} terms, c1' sum is known to be {sum_value}; error is {err(sum_value)}")
    return sum_value 

def c1p_const():
    return 6/e**2 - (log(2)**2 + log(2) + 3)/(2*e)

# the antiderivative of (y-n) (log y) (log(me) - log y) / y^2
def c1p_antideriv(n,m,a):
    ans = 0
    ans += n * log(m*e) * (log(a)+1)/a
    ans += log(m*e) * log(a)**2 / 2
    ans -= log(a)**3 / 3
    ans -= n * (log(a)**2 + 2*log(a) + 2) / a 
    return ans

# Computes the integral 1/e \int_e^infty {y} log (\lceil y/e \rceil/(y/e)) (log y) dy to high precision by evaluating the integral exactly
def c1p_integral(T):
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

        integral += c1p_antideriv(n,m,next_a)
        integral -= c1p_antideriv(n,m,a)

        if terminate:
            break

        a = next_a
        n = next_n
        m = next_m

    # estimate integrand by [0,e log y/y^3], so tail integral is in [0, e(2 log T+1)/4T^2]
    integral += mpi(0, e * (2*log(T) + 1) / (4*T**2))
    integral /= e

    print(f"On integrating up to {T}, c1' integral is known to be {integral}; error is {err(integral)}")
    
    return integral

def c1p(K,T):
    c1p_val = c1p_sum(K) + c1p_const() - c1p_integral(T)
    print(f"c1' = {c1p_val}; error is {err(c1p_val)}") 
    return c1p_val

# Evaluate c1'' = \sum_{k=1}^infty 1/k log( \lfloor k/e \rfloor / (k/e) ) to high precision, by summing K terms and rigorously bounding the error
def c1pp_easy(K):
    sum = mpi(0,0)
    for k in range(1, K + 1):
        term = log((e/k) * ceil(k/e)) / k
        sum += term
# bound summand crudely by e/k^2, then estimate sum by integral test
    sum += mpi(0, e/K) 
    print(f"After {K} terms, c1'' sum is known to be {sum}; error is {err(sum)}")
    return sum 

# same as previous sum, but now implementing the bounds from the Erdos-Turan inequality
def c1pp_advanced(K,N):
    sum = mpi(0,0)
    for k in range(1, K + 1):
        term = log((e/k) * ceil(k/e)) / k
        sum += term

    # deal with an error term from Taylor expansion of the logarithm
    sum -= mpi(0, e*e/(4*K*K))

    # deal with the term where {k/e} is approximated by 1/2
    sum += mpi(e/(2*(K+1)), e/(2*K))

    # deal with the first term in Erdos-Turan
    sum += mpi(-e/((N+1)*K), e/((N+1)*K))

    for n in range(1, N + 1):
        x = ((2*e)/(pi*n) + (2*e)/(N+1)) / (abs(sin(pi*n/e)) * (K+1)**2)
        sum += mpi(-x,x)

    print(f"After {K} terms and {N} frequencies, c1'' sum is known to be {sum}; error is {err(sum)}")
    return sum


def c1_easy(c0_val,K,T,K2):
    return c1p(K,T) + c0_val * c1pp_easy(K2) - c0_val * c0_val * e / 2

def c1_advanced(c0_val,K,T,K2,N):
    return c1p(K,T) + c0_val * c1pp_advanced(K2,N) - c0_val * c0_val * e / 2

c0_val = c0(100000, 100000)
print(f"c0 = {c0_val}; error is {err(c0_val)}")

c1_val = c1_advanced(c0_val, 1000000, 100000, 1000000, 100000) 
print(f"c1 = {c1_val}; error is {err(c1_val)}")
