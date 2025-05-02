import matplotlib.pyplot as plt
import math

# Code for generating Figure 5 in the new paper.

def is_prime(n):
    """Return True if n is a prime number, else False."""
    if n < 2:
        return False
    # Check divisibility from 2 up to sqrt(n)
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def f(x):
    return math.floor(1/x) * math.log(math.ceil(1/(math.e * x)) * math.e * x)

def contrib(N,p):
    if is_prime(p):
        return f(p/N)
    return 0

def upper(N):
    return sum(contrib(N,p) for p in range(math.ceil((N/math.e) / math.floor(math.sqrt(N/math.e))), N + 1) if is_prime(p))

def upper2(N):
    return sum(contrib(N,p) for p in range(math.ceil(N/math.e), N + 1) if is_prime(p))

def lower(N):
    return math.log(2*math.pi*N) + 1/(12*N)

base = range(80,5000)
LHS = [lower(N) for N in base]
RHS = [upper(N) for N in base]
LHS2 = [upper2(N) for N in base]

plt.figure(figsize=(8, 6))
plt.plot(base, LHS, label='LHS of (5.6)' )
plt.plot(base, LHS2, label='LHS of (5.7)' )
plt.plot(base, RHS, label='RHS' )
plt.title('Left and right-hand sides of (5.6), (5.7)')
plt.xlabel('$N$')
plt.legend()
plt.grid(True)
plt.show()

