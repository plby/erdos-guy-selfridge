import matplotlib.pyplot as plt
import math

# Code for plotting distance to next smooth number

def is_smooth(N):
    while N % 2 == 0:
        N /= 2
    while N % 3 == 0:
        N /= 3
    return N == 1

def next_smooth(x):
    N = math.ceil(x)
    while not is_smooth(N):
        N += 1
    return N

# See Table 2
def kappa(L):
    if L >= 341.34:
        return math.log(9/8)
    if L >= 40.5:
        return math.log(32/27)
    if L >= 4.5:
        return math.log(4/3)
    if L >= 4/3:
        return math.log(3/2)
    return math.log(2)


base = [i/60 for i in range(60,6000)]
kappa_plot = [math.log(next_smooth(x)/x) for x in base]
kappa_upper = [kappa(x) for x in base]


plt.figure(figsize=(8, 6))
plt.plot(base, kappa_plot, label='$\\log (\\lceil x \\rceil^{\\langle 2,3 \\rangle}/x)$' )
plt.plot(base, kappa_upper, label='$\\kappa_x$' )
plt.title('Distance to next 3-smooth number')
plt.xlabel('$x$')
plt.legend()
plt.grid(True)
plt.show()

