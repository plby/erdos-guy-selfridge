from collections import List
from math import ceil, exp, floor, log, sqrt
from time import perf_counter_ns

#
# utility program for calculating bounds on t(N) implied
# by lemma 5.1 of the manuscript
#

fn sieve(N: Int, mut p: List[Float64]) -> Int:
    print('Sieving primes... ', end='')
    var t0 = perf_counter_ns()
    var np = 0
    for i in range(2, N+1):
        if p[i] != 0.0:
            continue  # i is not prime
        j = i
        while j <= N:
            p[j] = 1 # mark j as composite
            j += i
        p[np] = Float64(i)
        np += 1
    print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
    return np


fn f(N: Float64, T: Float64, p: Float64) -> Float64:
    # equation (1.7) of the manuscript with
    #   x = p/N
    #   a = N/T
    return floor(N/p)*log(ceil(T/p)*p/T)


fn main():

    var Nmax = 900_000_000
    #var Nmax = 300_000
    var p: List[Float64] = List[Float64](length=Nmax+1, fill=0)

    var np = sieve(Nmax, p)

    for E in range(4, 9):
        for M in range(1, 10):
            var N = M*10**E
            print('Bisecting N=' + String(N) + '... ', end='')
            var t0 = perf_counter_ns()

            var L = 0.0
            for i in range(2, N+1):
                L += log(Float64(i))

            var T_lo = 2*N//7
            var T_hi = N//2
            while T_hi-T_lo > 1:
                var T = (T_hi+T_lo)//2

                var n = Float64(N)
                var t = Float64(T)

                var B = L - n*log(t)
                var p_min = t/floor(sqrt(t))
                
                var s = 0.0
                var ix = np-1
                while p[ix] > p_min and s <= B:
                    s += f(n, t, p[ix])
                    ix -= 1
                #print('***', T, B, s, ix)

                if s > B:
                    T_hi = T
                else:
                    T_lo = T

            print(String((perf_counter_ns()-t0)//1_000_000) + 'ms')
            print('Bound', N, T_lo)


