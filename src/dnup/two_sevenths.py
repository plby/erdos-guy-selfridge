# This program is supposed to verify that t_{2,3,5,7} >= N*2/7 for
# sufficiently large N.  It follows the recipe of Proposition 6.7.

from fractions import Fraction

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

# This is the asymptotic fraction we are targeting.
alpha = Fraction(2,7)

# This is our set \mathcal{D}, the 7-smooth numbers up to 28.
primes = [2,3,5,7]
D = [1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,25,27,28]

# This is the initial segment of the infinite sequence of non-negative
# real numbers a_\ell.
a = {
    2: Fraction(145686,1000000),
    3: Fraction(10165,100000),
    4: Fraction(721964,10000000),
    5: Fraction(592374,10000000),
    6: Fraction(512691,10000000),
    7: Fraction(411295,10000000),
    8: Fraction(351711,10000000),
    9: Fraction(30096,1000000),
    10: Fraction(240768,10000000),
    12: Fraction(361152,10000000),
    14: Fraction(257966,10000000),
    15: Fraction(103186,10000000),
    16: Fraction(902881,100000000),
    18: Fraction(15048,1000000),
    20: Fraction(120384,10000000),
    21: Fraction(515932,100000000),
    24: Fraction(128983,10000000),
    25: Fraction(361152,100000000),
    27: Fraction(642048,100000000),
    30: Fraction(802561,100000000),
    32: Fraction(334431,100000000),
    36: Fraction(86941,10000000),
    40: Fraction(497602,100000000),
    42: Fraction(362285,100000000),
    48: Fraction(494001,100000000),
    50: Fraction(331489,100000000),
    }

# As an infinite sequence, it continues with a[x]=a[x/2]/2 if x/2 >=
# repeat_from.  For example, a[52] = a[26]/2.  As such, we will have
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

# Check the identity (6.7)
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
# >= repeat_from, we will assume that min(1/d, alpha/ell)=alpha_ell.
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
