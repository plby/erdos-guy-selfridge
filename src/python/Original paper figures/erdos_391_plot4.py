import numpy as np
import matplotlib.pyplot as plt
import math

# code for generating Figure 1 in the original paper, covering the range N <= 80

def plot3():
    L = 15000
    base = [i/L for i in range(1, L+1)]
    fn = [math.floor(L/i) * math.log(math.e * i/L * math.ceil(L/(math.e*i))) / math.e for i in range(1,L+1)]
    compare = [0.3044019010 for i in range(1,L+1)]
    upper = [math.log(1 + math.e*i/L)/(math.e*i/L) for i in range(1,L+1)]
    basetrunc = [i/L for i in range(math.ceil(L/math.sqrt(2*math.e)), L+1)]
    trunccompare = [math.log(math.e/2)/math.e for _ in range(math.ceil(L/math.sqrt(2*math.e)), L+1)]
    plt.figure(figsize=(8, 6))
    plt.plot(base, fn, label='$\\frac{1}{e} f_e(x)$' )
    plt.plot(base, compare, label='$c_0$' )
#    plt.plot(basetrunc, trunccompare, label='$\\frac{1}{e} \\log(e/2)$' )
    plt.plot(base, upper, label='$\\frac{\\log(1+ex)}{ex}$' )
    plt.title('$\\frac{1}{e} f_e(x)$')
    plt.xlabel('$x$')
    plt.legend()
    plt.grid(True)
    plt.show()

plot3()
