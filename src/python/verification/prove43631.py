import math
from fractions import Fraction
from sympy import factorint, primerange

# For an explanation of what's going on here, see this comment by
# Terence Tao:
# https://terrytao.wordpress.com/2025/03/26/decomposing-a-factorial-into-large-factors/#comment-687665

# This program proves that t(N) >= T using a certificate obtained from
# linear programming.  It uses exact arithmetic.
N = 43631
T = int(math.ceil(N/3))
print(f"This program attempts to prove t({N}) >= {T}.")

# We will use the following bounds.  Note that they are made slightly
# less redundant by avoiding repetition.  For example, a_p will be 1
# for all primes beginning at p>=14549, which is recorded just as a
# single entry for 14549 instead of once for each prime.

# The bound is tight (though we don't prove that with a corresponding
# primal).

# All the numbers are given as fractions over 1257 = 3*419, so a_2 = 91/1257
# and a_3 = 144/1257 = 48/419.
denominator = 1257
lp_bounds = {
2: 91,
3: 144,
5: 211,
7: 255,
11: 315,
13: 336,
17: 371,
19: 385,
23: 411,
29: 441,
31: 450,
37: 474,
41: 487,
43: 493,
47: 505,
53: 520,
59: 536,
61: 538,
67: 551,
71: 559,
73: 562,
79: 572,
83: 578,
89: 588,
97: 600,
101: 605,
103: 607,
107: 613,
109: 615,
113: 620,
127: 635,
131: 638,
137: 644,
139: 647,
149: 656,
151: 657,
157: 663,
163: 667,
167: 669,
173: 675,
179: 679,
181: 681,
191: 687,
193: 690,
197: 692,
211: 702,
223: 707,
227: 710,
229: 711,
233: 714,
239: 719,
251: 725,
257: 728,
263: 729,
269: 731,
271: 734,
277: 737,
281: 739,
293: 744,
307: 749,
311: 752,
317: 755,
331: 760,
347: 767,
359: 770,
367: 773,
373: 777,
383: 781,
397: 783,
409: 787,
419: 791,
431: 795,
443: 798,
457: 802,
479: 807,
487: 811,
503: 816,
521: 820,
541: 825,
563: 830,
587: 835,
607: 840,
641: 846,
673: 851,
701: 858,
733: 864,
769: 872,
809: 878,
857: 886,
911: 893,
971: 902,
1039: 911,
1123: 921,
1213: 931,
1327: 942,
1459: 955,
1619: 969,
1823: 984,
2081: 1002,
2437: 1022,
2909: 1046,
3637: 1075,
4861: 1113,
7283: 1166,
14549: 1257,
}

# We will have a coefficient a_p for each prime p up to N, derived
# from the previous table.  All other (ie, larger) primes also get a
# coefficient of 1 by default.
primes = list(primerange(2, N+1))
a = dict()
latest_a = None
for p in primes:
    if p in lp_bounds:
        latest_a = Fraction(lp_bounds[p], denominator)
    a[p] = latest_a

# We will need a function that returns the valuation of n at p, that
# is the highest power of p that divides n.  This function assumes p
# is prime.
def valuation(n, p):
    v = 0
    t = n
    while (t % p) == 0:
        v += 1
        t //= p
    return v

# We will use a more efficient routine for the valuation of N!
def valuation_of_factorial(n, p):
    v = 0
    t = n
    while t > 0:
        t //= p
        v += t
    return v

# Now we compute the bound we get from this whole setup.
bound = 0
for p in primes:
    bound += a.get(p, Fraction(1,1)) * valuation_of_factorial(N, p)

assert bound < N, "The resulting bound is insufficiently strong."
print(f"Specifically, it proves a bound on how many factors >={T} one can find in {N}!.")
print(f"The upper bound achieved is {bound} ~ {float(bound)}, which is less than {N}.")

# Check that the a_p are positive.
for p in primes:
    assert 0 <= a[p]

# Check that the appropriate inequalities are true up to N
for j in range(T, N+1):
    lhs = 0
    for p in factorint(j):
        lhs += a.get(p, Fraction(1,1)) * valuation(j, p)
    assert lhs >= 1, f"LP inequality failed for j={j}: lhs={lhs} ~ {float(lhs)}."
print(f"The appropriate LP inequality was verified up to j <= N={N}.")

# Check that the a_p are sorted.
previous_a = Fraction(0, 1)
for p in primes:
    assert previous_a <= a[p], f"a_p are not sorted at p={p}."
    previous_a = a[p]
assert previous_a <= 1, f"a_{p} exceeds 1, but all large primes will have a_p=1."
print(f"The a_p are sorted in (weakly) increasing order.")
