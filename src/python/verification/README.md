# Some verification-related Python scripts

This directory contains two files:

* `simple_check.py` verifies the factorizations in
  `Data/factorizations`.  See the README in that directory for more
  information.

* `prove43631.py` checks the dual linear programming certificate for
  the claim that $t(N) < N/3$ for $N=43631$.  The following output is
  produced:

```
erdos-guy-selfridge$ time python3 src/python/verification/prove43631.py 
This program attempts to prove t(43631) < 14544.
Specifically, it proves a bound on how many factors >=14544 one can find in 43631!.
The upper bound achieved is 54844120/1257 ~ 43630.96260938743, which is less than 43631.
The appropriate LP inequality was verified up to j <= N=43631.
The a_p are sorted in (weakly) increasing order.

real	0m1.314s
user	0m1.251s
sys	0m0.044s
```
