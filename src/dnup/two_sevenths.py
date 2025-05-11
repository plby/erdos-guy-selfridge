# This program is supposed to verify that t_{2,3,5,7} >= N*2/7 for
# sufficiently large N.  It follows the recipe of Proposition 6.3.

from fractions import Fraction
import math

def valuation(N: int, P: int) -> int:
    """
    Compute the P-adic valuation of N.

    The P-adic valuation v_P(N) is the largest exponent k such that P**k divides N.
    """
    result = 0
    while (N%P) == 0:
        N //= P
        result += 1
    return result

# We hope to prove that everything works if N >= N0 here:
N0 = 9*10**6

# This is the asymptotic fraction we are targeting.
alpha = Fraction(2,7)

# This is our set \mathcal{D}, the 7-smooth numbers up to 28.
primes = [2,3,5,7]
D = [1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,25,27,28]

# This is the initial segment of the infinite sequence of non-negative
# real numbers a_\ell.
a = {
    2: Fraction(145624,1000000),
    3: Fraction(101606,1000000),
    4: Fraction(721655,10000000),
    5: Fraction(59212,1000000),
    6: Fraction(512472,10000000),
    7: Fraction(411119,10000000),
    8: Fraction(35156,1000000),
    9: Fraction(300831,10000000),
    10: Fraction(240665,10000000),
    12: Fraction(360998,10000000),
    14: Fraction(257856,10000000),
    15: Fraction(103142,10000000),
    16: Fraction(902494,100000000),
    18: Fraction(150416,10000000),
    20: Fraction(120333,10000000),
    21: Fraction(515711,100000000),
    24: Fraction(128928,10000000),
    25: Fraction(360654,100000000),
    27: Fraction(639716,100000000),
    28: Fraction(29255,10000000),
    30: Fraction(517733,100000000),
    32: Fraction(453017,100000000),
    36: Fraction(755028,100000000),
    40: Fraction(604022,100000000),
    42: Fraction(716716,1000000000),
    45: Fraction(482397,100000000),
    48: Fraction(351965,100000000),
    50: Fraction(181207,100000000),
    }

# As an infinite sequence, it continues with a[x]=a[x/2]/2 if x/2 >=
# repeat_from.  For example, a[54] = a[27]/2.  As such, we will have
# to sum infinite series below (manually).
repeat_from = 26

# Compute sigma_{d,D}
sigma = dict()
for d in D:
    sigma[d] = Fraction(1,1)
    # Compute P_+(d):
    Pplusd = 1
    for p in primes:
        if (d % p) == 0:
            Pplusd = p
    # Defining product for sigma_{d,D}
    for p in primes:
        if (p < Pplusd or p*d in D):
            sigma[d] *= Fraction(p-1, p)

# Check the identity (6.7).
test = 0
for d in D:
    test += sigma[d] * Fraction(1,d)
assert test == 1

# Budgets (6.8)
for p in primes:
    lhs = 0
    for ell in a:
        v = valuation(ell, p)
        if ell >= repeat_from:
            if p == 2:
                # This corresponds to an infinite series where the
                # valuation also goes up by 1 each time.
                v = 2*(v+1)
            else:
                # This corresponds to an infinite series where the
                # valuation is not increasing as we go.
                v = 2*v
        lhs += v * a[ell]
    rhs = 0
    for d in sigma:
        rhs += sigma[d] * Fraction(valuation(d, p), d)
    assert lhs <= rhs

# The following assertion is used in a technical way in verifying the
# next condition.  Specifically, when we verify the condition for ell
# >= repeat_from, we will assume that min(1/d, alpha/ell)=alpha/ell.
# So this checks that:
assert Fraction(1,max(D)) > alpha * Fraction(1/repeat_from)

# Earth-moving (6.9).
for ell in a:
    lhs = 0
    for ellp in a:
        # The conditions below are somewhat tricky in terms of summing
        # the correct infinite series.
        if ellp > ell:
            if ellp >= repeat_from:
                # Infinite series
                lhs += 2*a[ellp]
            else:
                lhs += a[ellp]
        if repeat_from <= ellp <= ell:
            # Infinite series, but corresponding to a[2*ellp], so this
            # is 2*a[2*ellp] = 2*(a[ellp]/2) = a[ellp].
            lhs += a[ellp]
    rhs = 0
    for d in sigma:
        here = sigma[d] * min(Fraction(1,d), alpha * Fraction(1,ell))
        rhs += here
    assert lhs >= rhs
print("Asymptotic conditions (6.8) and (6.9) verified:")
print(f"t(N) >= alpha N for alpha={alpha} and sufficiently large N.")
print("")

# Now we will try to estimate how "sufficiently large" N must be.

# We use the modified sequence a'_\ell = ceil(a_\ell N)/N for \ell <=
# 2^L*N, and 0 otherwise.  This means we are truncating the infinite
# series, and we will have to account for that accordingly below.
# This value was chosen by hand.
L = 7

# For ell >= repeat_from, figure out what the sum_{ell'>ell} a_{ell'}
# is like, namely for what value A is this sum <= A/ell?
A = 0
for ell in a:
    if ell < repeat_from:
        continue
    here = 0
    for ellp in a:
        # The conditions below are somewhat tricky in terms of summing
        # the correct infinite series.
        if ellp > ell:
            if ellp >= repeat_from:
                # Infinite series
                here += 2*a[ellp]
            else:
                assert False # Should have: ellp > ell >= repeat_from
        if repeat_from <= ellp <= ell:
            # Infinite series, but corresponding to a[2*ellp], so this
            # is 2*a[2*ellp] = 2*(a[ellp]/2) = a[ellp].
            here += a[ellp]
    if A < here*ell:
        A = here*ell

# Let's figure out how much slack there is in (6.3).  If ell >
# alpha*N, the identity is automatically true, so we may assume ell <=
# alpha*N.

# The LHS has two effects going on: a ceiling on a_{ell'} only makes
# it bigger, so that helps us.  However, the sum is truncated, so that
# makes it smaller.  We will estimate how much smaller it is using A.

# Since ell <= alpha*N, 2^L*N >= 2^L/alpha*ell.

# If sum_{ell'>ell} a_{ell'} <= A/ell, then sum_{ell'>2^L*N} >=
# sum_{ell'>(2^L/alpha)*ell} a_{ell'} <= (A*alpha/2^L)/ell.  So
# overall, the LHS goes down by at most this much.

# The RHS has an error because of the roughness conditions in the
# count.  For roughness, we use an estimate that the number of
# 7-smooth numbers up to x is x*8/35 + O_{\le}(53/35).  The equivalent
# estimates for 5-, 3- and 2-smooth are tighter.
smoothness_bound = Fraction(53,35)
for ell in a:
    lhs = 0
    for ellp in a:
        # The conditions below are somewhat tricky in terms of summing
        # the correct infinite series.
        if ellp > ell:
            if ellp >= repeat_from:
                # Infinite series
                lhs += 2*a[ellp]
            else:
                lhs += a[ellp]
        if repeat_from <= ellp <= ell:
            # Infinite series, but corresponding to a[2*ellp], so this
            # is 2*a[2*ellp] = 2*(a[ellp]/2) = a[ellp].
            lhs += a[ellp]
    rhs = 0
    smoothness_error_denom = 0
    for d in sigma:
        here = sigma[d] * min(Fraction(1,d), alpha * Fraction(1,ell))
        rhs += here
        smoothness_error_denom += smoothness_bound
    if ell >= repeat_from:
        lhs -= (A*alpha/2**L)/ell
    assert lhs >= rhs
    epsilon = lhs - rhs
    assert N0 > smoothness_error_denom/epsilon

# Now let's figure out how much slack there is in (6.2).  The RHS has
# the same smoothness issues as in the previous section.

# The LHS gets bigger because of all of the ceilings.  Each ceiling
# causes a local increase of at most 1/N.  Each term is at most
# valuation_p(ell) <= valuation_2(2^L*N) = (L + log_2(N)) if p=2, and
# at most the maximum valuation of the extreme small ell otherwise.

# Let's estimate the number of terms.
odd_ells = 0
for ell in a:
    if (ell % 2) == 1:
        odd_ells += 1
# Each ell in our LHS has one of these odd terms times a power of 2,
# so the number of terms is at most odd_ells*log_2(N).  Note that
# there's possibly (log_2(N)+1) terms for the odd ell 1 (this is a
# fencepost-counting issue), but there's much fewer for larger odd ell
# so they compensate.
assert odd_ells >= 2

# All in all, the total increase to the LHS is at most
# (L+log_2(N))*odd_ells*log_2(N)/N for p=2 and a simpler thing if
# p>2.  Either way, this is decreasing for N >= e^2.  Let's make sure
# that N is at least that big, even though of course it is..
assert N0 > math.exp(2)

for p in primes:
    lhs = 0
    for ell in a:
        v = valuation(ell, p)
        if ell >= repeat_from:
            if p == 2:
                # This corresponds to an infinite series where the
                # valuation also goes up by 1 each time.  (For
                # simplicity, we don't even truncate the sums!  This
                # only makes it worse for us.)
                v = 2*(v+1)
            else:
                # This corresponds to an infinite series where the
                # valuation is not increasing as we go.  (For
                # simplicity, we don't even truncate the sums!  This
                # only makes it worse for us.)
                v = 2*v
        lhs += v * a[ell]
    error_estimate = (L+math.log2(N0))*odd_ells*math.log2(N0)
    if p > 2:
        max_small_ell_valuation = max([valuation(small_ell, p) for small_ell in a])
        error_estimate = max_small_ell_valuation*odd_ells*math.log2(N0)
    lhs += error_estimate/N0

    smoothness_error_denom = 0
    rhs = 0
    for d in sigma:
        rhs += sigma[d] * Fraction(valuation(d, p), d)
        smoothness_error_denom += valuation(d, p)

    assert lhs <= rhs
    epsilon = rhs - lhs
    assert N0 > smoothness_error_denom/epsilon

print("Conditions (6.2) and (6.3) verified for the modified sequence:")
print(f"t(N) >= alpha N for alpha={alpha} and N >= {N0}.")
print("")
