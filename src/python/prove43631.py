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
lp_bounds = {
2: Fraction(723946, 10000000),
3: Fraction(1145585, 10000000),
5: Fraction(16786, 100000),
7: Fraction(202864, 1000000),
11: Fraction(2505967, 10000000),
13: Fraction(2673031, 10000000),
17: Fraction(2951472, 10000000),
19: Fraction(3062848, 10000000),
23: Fraction(326969, 1000000),
29: Fraction(3508353, 10000000),
31: Fraction(3579952, 10000000),
37: Fraction(3770883, 10000000),
41: Fraction(3874304, 10000000),
43: Fraction(3922037, 10000000),
47: Fraction(4017502, 10000000),
53: Fraction(4136834, 10000000),
59: Fraction(4264121, 10000000),
61: Fraction(4280032, 10000000),
67: Fraction(4383453, 10000000),
71: Fraction(4447096, 10000000),
73: Fraction(4470963, 10000000),
79: Fraction(4550517, 10000000),
83: Fraction(459825, 1000000),
89: Fraction(4677804, 10000000),
97: Fraction(477327, 1000000),
101: Fraction(4813047, 10000000),
103: Fraction(4828958, 10000000),
107: Fraction(4876691, 10000000),
109: Fraction(4892601, 10000000),
113: Fraction(4932379, 10000000),
127: Fraction(505171, 1000000),
131: Fraction(5075577, 10000000),
137: Fraction(5123309, 10000000),
139: Fraction(5147176, 10000000),
149: Fraction(5218775, 10000000),
151: Fraction(522673, 1000000),
157: Fraction(5274463, 10000000),
163: Fraction(5306285, 10000000),
167: Fraction(5322196, 10000000),
173: Fraction(5369928, 10000000),
179: Fraction(540175, 1000000),
181: Fraction(5417661, 10000000),
191: Fraction(5465394, 10000000),
193: Fraction(548926, 1000000),
197: Fraction(5505171, 10000000),
211: Fraction(5584726, 10000000),
223: Fraction(5624503, 10000000),
227: Fraction(5648369, 10000000),
229: Fraction(5656325, 10000000),
233: Fraction(5680191, 10000000),
239: Fraction(5719968, 10000000),
251: Fraction(5767701, 10000000),
257: Fraction(5791567, 10000000),
263: Fraction(5799523, 10000000),
269: Fraction(5815434, 10000000),
271: Fraction(58393, 100000),
277: Fraction(5863166, 10000000),
281: Fraction(5879077, 10000000),
293: Fraction(5918854, 10000000),
307: Fraction(5958632, 10000000),
311: Fraction(5982498, 10000000),
317: Fraction(6006364, 10000000),
331: Fraction(6046142, 10000000),
347: Fraction(610183, 1000000),
359: Fraction(6125696, 10000000),
367: Fraction(6149562, 10000000),
373: Fraction(6181384, 10000000),
383: Fraction(6213206, 10000000),
397: Fraction(6229117, 10000000),
409: Fraction(6260939, 10000000),
419: Fraction(6292761, 10000000),
431: Fraction(6324582, 10000000),
443: Fraction(6348449, 10000000),
457: Fraction(638027, 1000000),
479: Fraction(6420048, 10000000),
487: Fraction(645187, 1000000),
503: Fraction(6491647, 10000000),
521: Fraction(6523469, 10000000),
541: Fraction(6563246, 10000000),
563: Fraction(6603023, 10000000),
587: Fraction(66428, 100000),
607: Fraction(6682578, 10000000),
641: Fraction(673031, 1000000),
673: Fraction(6770088, 10000000),
701: Fraction(6825776, 10000000),
733: Fraction(6873508, 10000000),
769: Fraction(6937152, 10000000),
809: Fraction(6984885, 10000000),
857: Fraction(7048528, 10000000),
911: Fraction(7104216, 10000000),
971: Fraction(7175815, 10000000),
1039: Fraction(7247414, 10000000),
1123: Fraction(7326969, 10000000),
1213: Fraction(7406523, 10000000),
1327: Fraction(7494033, 10000000),
1459: Fraction(7597454, 10000000),
1619: Fraction(7708831, 10000000),
1823: Fraction(7828162, 10000000),
2081: Fraction(797136, 1000000),
2437: Fraction(8130469, 10000000),
2909: Fraction(83214, 100000),
3637: Fraction(8552108, 10000000),
4861: Fraction(8854415, 10000000),
7283: Fraction(9276054, 10000000),
14549: Fraction(1, 1),
}

# We will have a coefficient a_p for each prime p up to N, derived
# from the previous table.  All other (ie, larger) primes also get a
# coefficient of 1 by default.
primes = list(primerange(2, N+1))
a = dict()
latest_a = None
for p in primes:
    if p in lp_bounds:
        latest_a = lp_bounds[p]
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
print(f"The appropriate LP inequality was verified up to j<=N={N}.")

# Check that the a_p are sorted.
previous_a = Fraction(0, 1)
for p in primes:
    assert previous_a <= a[p], f"a_p are not sorted at p={p}."
    previous_a = a[p]
assert previous_a <= 1, f"a_{p} exceeds 1, but all large primes will have a_p=1."
print(f"The a_p are sorted in (weakly) increasing order.")
