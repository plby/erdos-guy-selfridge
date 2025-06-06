# Rearrangement-related utilities

This directory contains six programs related to Section
"6. Rearranging the standard factorization":

* `rearrange_2.cc` computes $t_2(N)$ very efficiently.  In particular,
  it can compute it for all $N$ up to $10^6$ in 1--2 seconds.  If you
  run the program, it just prints $t(N)$ for ever-increasing values of
  $N$.  As a result, we can clearly see the asymptotic behavior
  $t_2(N)/N \to 3/16 = 0.1875$:

```
erdos-guy-selfridge$ g++ -O2 -std=c++11 src/rearrange/rearrange_2.cc -o rearrange_2
erdos-guy-selfridge$ time ./rearrange_2 | head -n 1000000 | tail
999991 187501
999992 187502
999993 187502
999994 187502
999995 187502
999996 187503
999997 187503
999998 187503
999999 187503
1000000 187504

real	0m1.805s
user	0m1.329s
sys	0m2.059s
```

* `rearrange_2_gs.cc` also computes $t_2(N)$ very efficiently, but
  this time using the explicit formula given by Guy & Selfridge in
  their paper "Factorial Factorial n".  The output of this program
  agrees with the previous one for all $N \le 10^6$.

* `rearrange_2_gs.py` is a Python implementation of this same explicit
  formula given by Guy & Selfridge.  The output of this program agrees
  with the previous two for all $N \le 10^6$.

* `rearrange_3.cc` helps compute $t_{2,3}(N)$ using dynamic
  programming.  (The program `rearrange_lp.py` below is more general;
  this is just an extra check on the case $P=3$.)  The program must be
  run with two arguments, $N$ and $T$, and then the program prints $1$
  if $t_{2,3}(N) \ge T$ or $0$ if $t_{2,3}(N) < T$.  For example, we
  can check the last line of `Data/rearrange/up_to_3.txt` as follows:

```
erdos-guy-selfridge$ g++ -O2 -std=c++11 src/rearrange/rearrange_3.cc -o rearrange_3
erdos-guy-selfridge$ tail -n 1 Data/rearrange/up_to_3.txt
3 50000 12492
erdos-guy-selfridge$ time ./rearrange_3 50000 12492
1

real	0m12.459s
user	0m12.450s
sys	0m0.005s
erdos-guy-selfridge$ time ./rearrange_3 50000 12493
0

real	0m12.029s
user	0m12.023s
sys	0m0.004s
```

* `rearrange_lp.py` helps compute $t_{2,\ldots,P}(N)$.  The program is
  best run with three arguments: $P$, $N$, and $T$.  It outputs a
  linear program that can be solved with `lp_solve`.  If the program
  is feasible, then $t_{2,\ldots,P}(N) \ge T$; if the program is
  infeasible, then $t_{2,\ldots,P}(N) < T$.  For example, we can
  double-check the last line of `Data/rearrange/up_to_3.txt` and check
  the last line of `Data/rearrange/up_to_7.txt` as follows:

```
erdos-guy-selfridge$ tail -n 1 Data/rearrange/up_to_3.txt
3 50000 12492
erdos-guy-selfridge$ time python3 src/rearrange/rearrange_lp.py --ilp 3 50000 12492 | lp_solve | head -n 2

Value of objective function: 0

real	0m0.490s
user	0m0.451s
sys	0m0.049s

erdos-guy-selfridge$ time python3 src/rearrange/rearrange_lp.py --ilp 3 50000 12493 | lp_solve | head -n 2
This problem is infeasible

real	0m0.482s
user	0m0.431s
sys	0m0.061s

erdos-guy-selfridge$ tail -n 1 Data/rearrange/up_to_7.txt
7 50000 14530
erdos-guy-selfridge$ time python3 src/rearrange/rearrange_lp.py --ilp 7 50000 14530 | lp_solve | head -n 2

Value of objective function: 0

real	0m0.571s
user	0m0.545s
sys	0m0.041s

erdos-guy-selfridge$ time python3 src/rearrange/rearrange_lp.py --ilp 7 50000 14531 | lp_solve | head -n 2
This problem is infeasible

real	0m0.574s
user	0m0.518s
sys	0m0.074s
```

* `rearrange_lp_sol.py` helps verify the lower bound produced by the
  `lp_solve` output of a program produced by the previous program.  It
  checks all of the rearrangement conditions, so if it succeeds, the
  input is verified.  Optionally, it is also possible to output a file
  that can be checked by `simple_check.py`.  For example:

```
erdos-guy-selfridge$ time python3 src/rearrange/rearrange_lp.py --ilp 7 50000 14530 | lp_solve | python3 src/rearrange/rearrange_lp_sol.py 7 50000 14530
7 50000 14530

real	0m2.736s
user	0m2.001s
sys	0m0.099s

erdos-guy-selfridge$ time python3 src/rearrange/rearrange_lp.py --ilp 7 50000 14530 | lp_solve | python3 src/rearrange/rearrange_lp_sol.py --simple_check 7 50000 14530 | python3 src/python/verification/simple_check.py
50000 50000 14530

real	0m6.608s
user	0m5.149s
sys	0m0.117s
```
