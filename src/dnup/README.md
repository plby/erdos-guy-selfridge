# Verifications for rearrangement-related computer-assisted proofs

This directory contains three scripts, and one associated file, to
verify several statements from the paper.  The techniques involved are
from Section "6. Rearranging the standard factorization".

In all cases, to perform the verification, you need to simply run the
file using `python3`.

* `verify.py` is a rearrangement-based proof that $t(N) \ge N/3$ for
  sufficiently large $N$.  No attempt is made to determine an
  effective bound on $N$ for which the statement holds.  The file
  `data.py` contained some associated data, in the form of (rounded)
  (dual) linear programming weights.

* `one_fourth.py` proves that $t_{2,3} < N/4$ for sufficiently large
  $N$.  An effective bound is obtained: $N \ge 1328148$.

* `two_sevenths.py` proces that $t_{2,3,5,7} \ge 2N/7$ for
  sufficiently large $N$.  The asymptotic and effective bounds are
  checked separately; the effective bound verified is $N \ge 9000000$.

The programs produce the following output:

```
erdos-guy-selfridge$ time python3 src/dnup/verify.py 
Extra powers of 2 available: ~ 0.00014473771961589284
Slack available in multiples: ~ 1.414436900160613e-06
Slack available in budget:  1776011/1215000000000 ~ 1.4617374485596707e-06
All asserts passed!  Down-reorganize-up algorithm succeeded.

real	0m4.700s
user	0m4.639s
sys	0m0.060s

erdos-guy-selfridge$ time python3 src/dnup/one_fourth.py 
We have proven that, for sufficiently large N, it is NOT possible to split N! into N factors,
each of which is >=N/4, by only moving factors of 2 and 3.

Specifically, it is not possible to get asymptotically more than b*N factors, where
b = 4457832185537/4458050224128 ~ 0.9999510910420389.

Note that b = 1 - epsilon, where epsilon ~ 1/20446.

The 'sufficiently large' can be quantified.  We get an inequality 1 <= b + O_\le(c / N) for b as above and
c = 1559/24 ~ 65, from which it follows that it holds for N >= c/(1-b):

N >= 1328148.

real	0m0.190s
user	0m0.183s
sys	0m0.005s

erdos-guy-selfridge$ time python3 src/dnup/two_sevenths.py 
Asymptotic conditions (6.8) and (6.9) verified:
t(N) >= alpha N for alpha=2/7 and sufficiently large N.

Conditions (6.2) and (6.3) verified for the modified sequence:
t(N) >= alpha N for alpha=2/7 and N >= 9000000.


real	0m0.065s
user	0m0.062s
sys	0m0.000s
```