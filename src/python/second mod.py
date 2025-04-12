import numpy as np
import matplotlib.pyplot as plt
import math

# code for generating a second modification of the c0 constant. 

K = 342



def integral():
    sum = 0
    for k in range(2,3*K+1):
        sum += (k-1) * (1/(k-1)*math.log(1/(k-1)) - 1/(k)*math.log(1/(k)) + (math.log(3*math.ceil(k/3))-1)*(1/(k-1)-1/(k)))
    return sum

def total_variation():
    sum = 0
    for k in range(2, 3*K+1):
        a = (k-2)*math.log((3/(k-1))*math.ceil((k-1)/3))
        b =  (k-1)*math.log((3/(k-1))*math.ceil(k/3))
        c = (k-1)*math.log((3/k)*math.ceil(k/3))
#        print(a,b,c)
        sum += abs(a-b) + abs(b-c)
    return sum



def plot3():
    D = 5000  # discretization parameter
    L = D*3*K

    base = [i/L for i in range(D, L+1)]
    fn = [math.floor(L/i) * math.log(3*i/L * math.ceil(L/(3*i))) for i in range(D,L+1)]
    
    sum = 0
    tv = 0
    prev = 0
    for s in fn:
        sum += s
        tv += abs(s-prev)
        prev = s
    tv += s
    c = sum / L
    
    print(f"Numerical integral = {c}")
    print(f"Exact integral = {integral()}")

    print(f"Numerical total variation = {tv}")
    print(f"Total variation = {total_variation()}")

    fullbase = [i/L for i in range(1,L+1)]
    compare = [c for _ in range(1,L+1)]

    plt.figure(figsize=(8, 6))
    plt.plot(base, fn, label='$\\lfloor \\frac{1}{x} \\rfloor \\log(3x \\lceil \\frac{1}{3x} \\rceil)$' )
    plt.plot(fullbase, compare, label='$c_1$' )
    plt.title('$\\lfloor \\frac{1}{x} \\rfloor \\log(3x \\lceil \\frac{1}{3x} \\rceil)$')
    plt.xlabel('$x$')
    plt.legend()
    plt.xlim(0,1)
    plt.grid(True)
    plt.show()

plot3()
