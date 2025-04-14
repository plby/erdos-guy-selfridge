from fractions import Fraction
import logging
from sympy import primerange

logging.basicConfig(level=logging.INFO)
#logging.basicConfig(level=logging.DEBUG)

# Parameters for the algorithm; note that we will later use that Kup
# is specifically 2**V.
V = 11
Kdn = 2**V
Kup = 2**V
assert Kup >= Kdn/3
kdn_range = list(range(1, Kdn+1))
kup_range = list(range(1, Kup+1))

# Load some external data, indicating how much to multiply by
from data import multiples
assert sorted(multiples.keys()) == kup_range

# Collect relevant primes
primes = list(primerange(2, max(Kdn, Kup)+1))
logging.debug(f"{primes=}")

## DOWN

# For the down part, we start with our list of integers L_1 = (0,N].
# We start with the prime 2 and split it into the even and odd
# elements.  We divide each even element by 2, and as a result we
# collect N/2 powers of 2.  We are left with the odd numbers in (0,N]
# and also, separately, all numbers in (0,N/2].  We will call the
# first list (which we will continue to change) L_1 and the second one
# L_2.  We continue by splitting off L_2 into L_4 and so forth, except
# we never go to an index beyond Kdn.  In the end, we will have: L_1
# is half of (0,N], L_2 is half of (0,N/2], ... and finally L_V is all
# of (0,N/2**V].  We also collect a bunch of powers of 2.

# Then we consider the prime 3.  The first step would be to take L_1
# and remove the multiples of 3.  L_1 is then left with 1/3 of the
# numbers in (0,N] while L_3 becomes half of the numbers in (0,N/3];
# it is only half because recall that L_1 was already missing the even
# numbers.  Also, we only get N/6 multiples of 3 as a result.  And so
# on for other uses of 3 and larger primes.

# We will collect the powers of primes in a "budget", where something
# like 1/2 will mean 1/2 of N.
budget = {p: Fraction(0, 1) for p in primes}

# We will also track the fraction of each list that we have.  The list
# L_k will contain fraction[k] of (0,N/k].
fraction = dict()
fraction[1] = Fraction(1, 1)

# Main loop
for p in primes:
    for k in sorted(fraction.keys()):
        # Keep trying to remove the multiples of p in k, kp, ...
        t = k
        while t*p <= Kdn:
            # We add to our budget
            budget[p] += Fraction(1,t) * fraction[t]/p
            
            # L_{tp} acquires the same fraction as L_t
            fraction[t*p] = fraction[t]

            # But now remove a fraction 1/p from L_t
            fraction[t] *= 1 - Fraction(1, p)

            # And keep going
            t *= p

assert sorted(fraction.keys()) == kdn_range, "Expected to collect lists L_1, ..., L_{Kdn}"

logging.debug(f"{budget=}")
logging.debug(f"{fraction=}")

## RE-ORGANIZATION

# Now we have a bunch of lists L_k (as well as a budget of different
# prime factors, which we will not touch here).  We will now
# re-organize these lists, grouping them by the least multiple needed
# for them to get up to 1/3.

# For example, L_1 consists of some fraction of (0,N].  We will
# re-group this list into sublists S_1 = (N/3, N], then S_2 = (N/6,
# N/3], then S_3 = (N/9, N/6], etc.  This is because the first sublist
# S_1 is already at least N/3; the second sublist S_2 needs to be
# multiplied by 2; the third S_3, by 3; and so forth.  We stop this
# process at Kup, so the final sublist like this is S_{Kup} =
# ((N/3)/Kup, (N/3)/(Kup-1)] and we have one extra special sublist
# S_{Kup+1} = (0, (N/3)/Kup].

# When we next consider L_2, which is some fraction of (0,N/2], we
# will take the sublist (N/3,N/2] and add it to the sublist S_1.  That
# is, we no longer track which L_k something came from.  We simply
# track what fraction of (0,N] we are accumulating.

# We will collect information in two ways.  If we have a sublist like
# (N/3, N/2], this will be counted as a partial interval and its
# fraction times length recorded in part[1].  If we have a sublist
# like (0, N/3], then this is a whole interval and its fraction will
# be recorded in whole[2].
part  = {k: Fraction(0, 1) for k in kup_range}
whole = {k: Fraction(0, 1) for k in kup_range}

for k in kdn_range:
    multiple = (k // 3)+1
    if (k % 3) == 0:
        # Just a whole interval.
        whole[multiple] += fraction[k]
    elif (k % 3) == 1:
        # One part (N/(k+2), N/k]
        length = Fraction(1, k) - Fraction(1, k+2)
        part[multiple] += length * fraction[k]
        # and one whole.
        whole[multiple+1] += fraction[k]
    else:
        # One part (N/(k+1), N/k]
        length = Fraction(1, k) - Fraction(1, k+1)
        part[multiple] += length * fraction[k]
        # and one whole.
        whole[multiple+1] += fraction[k]

# Now we will add up the total fractions in each sublist, so total[k]
# denote the total fraction of (0, N] that ended up in S_k.  The
# variable whole_sum keeps track of how many whole intervals we have
# accumulated so far.
total = {k: Fraction(0, 1) for k in kup_range}
whole_sum = 0
for k in kup_range:
    whole_sum += whole[k]
    if k == 1:
        # Just avoiding a division by 0 here, which shouldn't come up
        # because the whole interval corresponding to k=1 would have
        # been (N/1, N/0].
        assert whole[k] == 0
        length = 0
    else:
        # The interval in question is ((N/3)/k), (N/3)/(k-1))].
        length = Fraction(1, 3*(k-1)) - Fraction(1, 3*k)

    total[k] = part[k] + length * whole_sum

# As we mentioned, there is one special sublist at the end corresponding to
# (0, (N/3)/Kup)
length = Fraction(1, 3*Kup)
special = length * whole_sum
logging.debug(f"{special=}")

# Now let's compute a cdf for the totals
total_cdf = {k: Fraction(0, 1) for k in kup_range}
total_cdf[0] = Fraction(0, 1)
for k in kup_range:
    total_cdf[k] = total_cdf[k-1] + total[k]
    logging.debug(f"total_cdf[{k}]={total_cdf[k]}")
assert total_cdf[Kup] + special == 1

# And compute the cdf for the multiples
multiples_cdf = {k: Fraction(0, 1) for k in kup_range}
multiples_cdf[0] = Fraction(0, 1)
for k in kup_range:
    multiples_cdf[k] = multiples_cdf[k-1] + multiples[k]
    logging.debug(f"multiples_cdf[{k}]={multiples_cdf[k]}")

# Fix up the last entry of multiples so that the total sum is correct
multiples[Kup] = total_cdf[Kup] - multiples_cdf[Kup-1]
multiples_cdf[Kup] = multiples_cdf[Kup-1] + multiples[Kup]

## UP

# Recall that the special interval is (0, (N/3)/Kup].  Everything
# there needs to be multiplied by strictly more than Kup = 2**V.  We
# will use powers of 2 to do it all.  So, half of the interval we will
# multiply by V+1, a quarter by V+2, and so forth.  The sum of
# (V+i)/2^i for i from 1 to infinity is V+2.  That is, on average we
# will multiply each element by V+2 powers of 2.  So that's our
# special spending.  We will make sure to account for it in our
# budget.
special_spent = special * (V+2)

# Now we verify that each sublist gets multiplied up by (at least) the
# right amount.  multiples[1] says what fraction of the original (0,N]
# list gets multiplied by 1, multiples[2] by 2, and so forth.  So we
# need to check that the partial sums of multiples are not greater
# than the partial sums of total, both of which we already computed.
epsilon1 = Fraction(1,1)
for k in kup_range:
    assert multiples_cdf[k] <= total_cdf[k]
    if k != Kup and epsilon1 > total_cdf[k] - multiples_cdf[k]:
        epsilon1 = total_cdf[k] - multiples_cdf[k]
assert multiples_cdf[Kup] == total_cdf[Kup]

# Now we verify that we didn't overspend our budget
epsilon2 = Fraction(1, 1)
for p in primes:
    spent = Fraction(0, 1)
    if p == 2:
        # Account for the special spending
        spent += special_spent
    for k in kup_range:
        # Compute the number of factors of p in k
        v = 0
        t = k
        while (t % p) == 0:
            v += 1
            t /= p
        # Add up our spending
        spent += v * multiples[k]
    assert spent <= budget[p]
    if p == 2:
        print("Extra powers of 2 available:", "~", float(budget[p] - spent))
    # Compute how much slack we have
    if epsilon2 > budget[p] - spent:
        epsilon2 = budget[p] - spent
        
print("Slack available in multiples:", "~", float(epsilon1))
print("Slack available in budget: ", epsilon2, "~", float(epsilon2))
print("All asserts passed!  Down-reorganize-up algorithm succeeded.")
