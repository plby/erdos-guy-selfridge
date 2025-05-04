from fractions import Fraction
import logging
import math

logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

# Parameters for the verification
MAX_TWO = 19
MAX_THREE = 12
M = 2**2 * 3**9 # 78732

# Generate some 3-smooth numbers, including all of them less than 4*M
assert 2**MAX_TWO > 4*M
assert 3**MAX_THREE > 4*M

smooth = []
for two in range(MAX_TWO+1):
    for three in range(MAX_THREE+1):
        x = 2**two * 3**three
        smooth.append(x)
smooth.sort()
logging.debug(f"{smooth=}")

# One-third of numbers are coprime to 2*3:
coprime_fraction = Fraction(1, 3)

# Figure out what fraction of factors require being multiplied by
# exactly x after all the powers of 2 and 3 have been removed.  This
# means the factor is in N/4 * [1/x, 1_xprev] where x_prev is the
# previous smooth number.  Recaling [0,N] to [0,1] from now on, this
# is the interval [1/(4*x), 1/(4*x_prev)].  We store the fraction of
# numbers that require multiplying by x in a dictionary dist (standing
# for distribution).
x_prev = None
dist = dict()
for x in smooth:
    if x > M:
        break
    dist[x] = Fraction(0, 1)
    # Suppose BEFORE dividing out 2s and 3s, the 3-smooth part of the
    # factor was y.  That means that BEFORE dividing out the 2s and
    # 3s, the factor was originally in the interval [y/(4*x),
    # y/(4*x_prev)].  In particular, if y >= 4*x, then this is not
    # possible.
    for y in smooth:
        if y >= 4*x:
            break
        # Now we compute the interval [y/(4*x), y/(4*x_prev)]
        # intersected with [0,1]
        lo = Fraction(y, 4*x)
        if x_prev is None or Fraction(y, 4*x_prev) > 1:
            hi = Fraction(1, 1)
        else:
            hi = Fraction(y, 4*x_prev)
        assert lo < hi
        assert hi <= 1

        # Add a 1/y fraction of the coprime fraction of our interval
        # to the distribution.  Note that y/(4*x) divided by y is
        # 1/(4*x), but the hi endpoint is trickier due to the
        # interaction with 1.
        dist[x] += coprime_fraction * (hi - lo) * Fraction(1, y)

    # Set x_prev for the next time around
    x_prev = x

# Let a_x be the fraction of numbers that are multiplied by x.  We
# will sum up a bunch of inequalities of the form sum_x lhs_x a_x <=
# rhs for different constants lhs_x and rhs.  The hope is to obtain
# something like x_1 + x_2 + ... < 1.  We will also multiply the
# inequalities by specially-chosen non-negative coefficients.
lhs = {x: Fraction(0, 1) for x in smooth}
rhs = Fraction(0, 1)

# The first two inequalities will be the prime budgets.  These are
# inequalities where the lhs is of the form sum_x v_p(x) a_x.

# All coefficients will have denominator 32.
denominator = 32

# The budget for 2 is 1.  The coefficient of this inequality is 2/32 =
# 1/16.
coefficient = Fraction(2, denominator)
rhs += coefficient * Fraction(1, 1)
for two in range(MAX_TWO+1):
    for three in range(MAX_THREE+1):
        x = 2**two * 3**three
        lhs[x] += coefficient * two

# The budget for 3 is 1/2.  The coefficient is 3/32.
coefficient = Fraction(3, denominator)
rhs += coefficient * Fraction(1, 2)
for two in range(MAX_TWO+1):
    for three in range(MAX_THREE+1):
        x = 2**two * 3**three
        lhs[x] += coefficient * three

# We don't have to check this here, but note that these coefficients
# have been chosen so that every x > M has a coefficient of at least 1
# already.
for x in smooth:
    if x > M:
        assert lhs[x] >= Fraction(1, 1)

# Next, we will add up some earth-moving constraints.  The
# earth-moving constraint for x says that a_1 + ... + a_x is bounded
# above by the proportion of factors that need to be multiplied by x
# or smaller.
for two in range(MAX_TWO+1):
    for three in range(MAX_THREE+1):
        x = 2**two * 3**three
        # The coefficient is 2/32 if x=1, 1/32 if x is nice, and 0
        # otherwise.  Nice means that x<=M and the power of 2 is at
        # most 2.
        coefficient = 0
        if two <= 2:
            coefficient = Fraction(1, denominator)
        if x == 1:
            coefficient = Fraction(2, denominator)
        if x > M:
            coefficient = 0
        if coefficient == 0:
            continue
        # Let's figure out the amount of earth available to move.
        # This is the cdf of our distribution up to x.
        cdf = 0
        for y in smooth:
            if y <= x:
                cdf += dist[y]
        # So that's the rhs of the inequality.
        rhs += coefficient * cdf

        # The lhs is just the partial sum a_1 + ... a_x:
        for y in smooth:
            if y <= x:
                lhs[y] += coefficient * Fraction(1, 1)

# Earlier we said that the goal was the lhs is of the form a_1 + ...,
# but actually all we need is that each coefficient is AT LEAST 1.  So
# we check that now.
for x in smooth:
    logging.debug(f"{x=} {float(lhs[x])=} {lhs[x]=}")
    assert lhs[x] >= Fraction(1, 1)

# And we want the rhs to be strictly less than 1, so our last step is
# to check that here.
logging.debug(f"{float(rhs)=} {rhs=}")
assert rhs < Fraction(1, 1)

# That completes the proof, but let's also compute the bound in
# another (faster) way and check that it's equal.  Along the way we
# also compute a bound on N when this argument takes effect.
def sumd(k):
    result = 0
    for m in smooth:
        if m < 4*k:
            result += (Fraction(1,m)-Fraction(1,4*k)) * Fraction(1,3)
    return result
def sume(k):
    result = 0
    for m in smooth:
        if m < 4*k:
            result += Fraction(4,3)
    return result

c = 0
b = 0
b += Fraction(2, denominator) * Fraction(1, 1) # 2-budget
b += Fraction(3, denominator) * Fraction(1, 2) # 3-budget
for two in range(MAX_TWO+1):
    for three in range(MAX_THREE+1):
        x = 2**two * 3**three
        # The coefficient is 2/32 if x=1, 1/32 if x is nice, and 0
        # otherwise.  Nice means that x<=M and the power of 2 is at
        # most 2.
        coefficient = 0
        if two <= 2:
            coefficient = Fraction(1, denominator)
        if x == 1:
            coefficient = Fraction(2, denominator)
        if x > M:
            coefficient = 0
        if coefficient == 0:
            continue
        b += coefficient * sumd(x)
        c += coefficient * sume(x)
assert rhs == b

print("We have proven that, for sufficiently large N, it is NOT possible to split N! into N factors,")
print("each of which is >=N/4, by only moving factors of 2 and 3.")
print("")
print("Specifically, it is not possible to get asymptotically more than b*N factors, where")
print(f"b = {b} ~ {float(b)}.")
print("")

print(f"Note that b = 1 - epsilon, where epsilon ~ 1/{math.floor(1/(1-b))}.")
print("")

print("The 'sufficiently large' can be quantified.  We get an inequality 1 <= b + O_\le(c / N) for b as above and")
print(f"c = {c} ~ {math.ceil(c)}, from which it follows that it holds for N >= c/(1-b):")
print("")

N0 = math.ceil(c/(1-b))
print(f"N >= {N0}.")
