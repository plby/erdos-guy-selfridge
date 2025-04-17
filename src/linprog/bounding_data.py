from find_crossover import log_random
from facfac import sieve, lpfac, greedy
import math
import matplotlib.pyplot as plt
import numpy as np

LOW = 10**3 # low end of intitial range
HIGH = 10**5 # high end of initial range
SAMPLE_COUNT = 300 # number of samples to take

def get_bounds(N, iteration_limit=5):
    # get initial upper and lower bounds on T(N)
    Tlb = math.ceil(2*N//7)
    if N == 56:
        Tlb -= 1
    Tub = math.ceil(N/math.exp(1))

    iterations = 0
    while Tub-Tlb > 1 and iterations < iteration_limit:
        T = (Tlb+Tub)//2
        print(T, Tlb, Tub)

        prob = sieve(N, T)

        fact, upper_bound = lpfac(prob)
        if N - upper_bound > 0.001: # check to make sure the upper bound for N is a decent bit below N
            Tub = T
            print("Shifting upper bound")
        else:
            fNT = greedy(prob, fact) # attempt to factorize
            nfactors = sum(fNT.f.values())
            if nfactors >= N:
                Tlb = T
                print("Shifting lower bound")
            else: # cannot conclusively establish a lower bound or an upper bound on this value
                break
        
        iterations += 1
    print(T, Tlb, Tub)

    return Tlb, Tub

samples = [log_random(LOW, HIGH) for i in range(SAMPLE_COUNT)] # sample points from log uniform distribution
samples.append(LOW) # add lower and upper points to the distribution
samples.append(HIGH)
samples.sort()

ITERATIONS = 10
bounds = [get_bounds(sample, iteration_limit=ITERATIONS) for sample in samples]

print(samples)
print(bounds)

samples = np.array(samples)

bounds = [(bounds[i][0]/samples[i], bounds[i][1]/samples[i]) for i in range(len(samples))]
bounds = np.array(bounds)

c0 = 0.30441901087
predicted_x = np.linspace(samples[0], samples[-1], num=1000)
predicted_y = 1/np.exp(1) - c0 / np.log(predicted_x)

# Plot the t(N)/N data along with the theoretical upper bound
plt.figure(figsize=(15, 5))
plt.scatter(samples, bounds[:, 1], label='Upper Bound')
plt.scatter(samples, bounds[:, 0], label='Lower Bound')
plt.plot(predicted_x, predicted_y, label=r'$y = \frac{1}{e} - \frac{c}{\ln(x)}$')
plt.xscale('log')
plt.xlabel('$N$ (log scale)')
plt.ylabel('$t(N)/N$')
plt.title('$t(N)/N$ Bounds')
plt.legend()
plt.tight_layout()
plt.savefig("t_large_n_bounds.png")
plt.show()

# Plot the difference between the theoretical upper bound and the actual data
new_bounds = np.copy(bounds)
new_bounds[:, 0] = (1/np.exp(1) - bounds[:, 1]) * np.log(samples) - c0
new_bounds[:, 1] = (1/np.exp(1) - bounds[:, 0]) * np.log(samples) - c0

print(new_bounds)

plt.figure(figsize=(15, 5))
plt.scatter(samples, new_bounds[:, 1], label='Upper Bound')
plt.scatter(samples, new_bounds[:, 0], label='Lower Bound')
plt.xscale('log')
plt.xlabel('$N$ (log scale)')
plt.ylabel('Estimated $o(1)$ Term')
plt.title('$o(1)$ Bounds')
plt.legend()
plt.tight_layout()
plt.savefig("o(1)_large_n_bounds.png")
plt.show()