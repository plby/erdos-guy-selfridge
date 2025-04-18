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

base = [i/60 for i in range(60,6000)]
base45 = [i/60 for i in range(270, 6000)]
base405 = [i/60 for i in range(2430, 6000)]
kappa = [math.log(next_smooth(x)/x) for x in base]
kappa45 = [math.log(4/3) for _ in base45]
kappa405 = [math.log(32/27) for _ in base405]


plt.figure(figsize=(8, 6))
plt.plot(base, kappa, label='$\\log (\\lceil x \\rceil^{\\langle 2,3 \\rangle}/x)$' )
plt.plot(base45, kappa45, label='$\\kappa_{4.5}$' )
plt.plot(base405, kappa405, label='$\\kappa_{40.5}$' )
plt.title('Distance to next 3-smooth number')
plt.xlabel('$x$')
plt.legend()
plt.grid(True)
plt.show()

