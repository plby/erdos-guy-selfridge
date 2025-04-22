import numpy as np
import matplotlib.pyplot as plt
import math

# code for generating a modification of the c0 constant.  Unfortunately, its average is too large for applications

def is_smooth(n):
    m = n
    while m % 2 == 0:
        m //= 2
    while m % 3 == 0:
        m //= 3
    if m == 1:
        return True
    return False

def roundup(n):
    m = math.ceil(n)
    while not is_smooth(m):
        m += 1
    return m

def plot3():
    L = 5000
    base = [i/L for i in range(1, L+1)]
    fn = [math.floor(L/i) * math.log(math.e * i/L * roundup(L/(math.e*i))) / math.e for i in range(1,L+1)]
    fn2 = [math.floor(L/i) * math.log(math.e * i/L * math.ceil(L/(math.e*i))) / math.e for i in range(1,L+1)]
    sum = 0
    for i in range(1,L+1):
        sum += fn[i-1]
    c = sum / L
    
    compare = [c for _ in range(1,L+1)]
    plt.figure(figsize=(8, 6))
    plt.plot(base, fn, label='$\\frac{1}{e} \\lfloor \\frac{1}{x} \\rfloor \\log(ex \\lceil \\frac{1}{ex} \\rceil^{2,3})$' )
    plt.plot(base, fn2, label='$\\frac{1}{e} \\lfloor \\frac{1}{x} \\rfloor \\log(ex \\lceil \\frac{1}{ex} \\rceil)$' )
    plt.plot(base, compare, label='$c_0$' )
    plt.title('$\\frac{1}{e} \\lfloor \\frac{1}{x} \\rfloor \\log(ex \\lceil \\frac{1}{ex} \\rceil)$')
    plt.xlabel('$x$')
    plt.ylim(0,2)
    plt.legend()
    plt.grid(True)
    plt.show()

plot3()
