# Factorizations supporting $t(N) \ge N/3$ for $N\in[43632, 80973]$.

One of the objectives of the paper is to prove the conjecture that
$t(N) \ge N/3$ for $N\ge 43632$.  For different ranges of $N$, this is
achieved in different ways.  This directory proves the conjecture for
$43632 \le N \le 80973$ by explicitly exhibiting factorizations.

The program `src/python/verification/simple_check.py` verifies that a
putative factorization of $N!$ into $F$ factors, each of which is at
least $M$, is correct.

The file format is a bunch of whitespace-separated integers. The first
three integers are $N$, $F$, and $M$ in that order.  Then there should
be $F$ integers which are the factors themselves, in any order, each
of which should be at least $M$.

The file `20-20-6.txt` in this directory is an example file for $20!$
factored as $20$ terms, each of which is at least $6$.  It contains
the following information, though with more newlines:

```
20 20 6
6 6 6 6 6 6 6 6
7 7
8 8
10 10 10 10
11 13 17 19
```

Similarly, the file `43632-43632-14545.txt` contains a factorization
of $43632!$ into $43632$ terms, each of which is at least $14545$.
Note that this proves that $t(N) \ge N/3$ (and in fact, $t(N)\ge
N/3+1$) for $N=43632$.

The remaining files, `43636-43636-14546.txt.gz` through
`80392-80392-26991.txt.gz`, prove similar results for larger values of
$N$.  Combined, they demonstrate that $t(N)\ge N/3$ for the range
$N\in[43632, 80973]$.  They use the fact that $t(N+1)\ge t(N)$ to
minimize the number of files involved.

The remaining files are compressed (note the `.gz` filenames).  The
script `simple_check.py` handles this transparently.  That is, you can
verify the files using commands like

```
erdos-guy-selfridge$ time python3 src/python/verification/simple_check.py Data/factorizations/43632-43632-14545.txt 
43632 43632 14545

real	0m0.982s
user	0m0.973s
sys	0m0.008s

erdos-guy-selfridge$ time python3 src/python/verification/simple_check.py Data/factorizations/80392-80392-26991.txt.gz 
80392 80392 26991

real	0m3.515s
user	0m3.467s
sys	0m0.028s
```

A couple of notes:

* The program `simple_check.py` makes no attempt to be efficient.

* The factorizations themselves were obtained via the integer
  programming technique, using the programs
  `src/linprog/generate_lp_format.cc` and `lp_solve`.

* An attempt was made to minimize the number of files involved in the
  verification in this directory, but it was not done perfectly.
