This directory contains a fast, memory efficient implemention of a variant of the greedy algorithm for constructing factorizations of $$N!$$ consisting of integers greater than $$t \ge \frac{N}{3}$$, with the goal of obtaining at least $$N$$ factors (thereby proving $$t(N) \ge t$$).  It is specifically optimized for large values of $$N$$, say $$N \ge 10^6$$; for smaller values of $$N$$ the bounds on $$t$$ it produces will typically not be as good as those given by the standard greedy algorithm (an implementation of which is also included), but it is much faster when $$N$$ is large (it is designed to handle $$N$$ up to $$10^12$$), and it uses only $$\tilde O(N^{\frac{5}{8}})$$ bits of memory (the 5/8 can be replaced by any exponent strictly above 1/2).

An implementation of the standard greedy algorithm is also included, but it's range is restricted to $$N\le 10^9$$ to keep the memory requirement low (under 1GB).

This implemention depends on Kim Walisch's [primesieve](https://github.com/kimwalisch/primesieve) and [primecount](https://github.com/kimwalisch/primecount) libraries for enumerating and counting primes, which must be installed before you compile the code in this directory.

These two libraries are used to efficiently handle factors divisible by primes in the range $$[\sqrt{t},N]$$. In this range the greedy algorithm can always use the optimal cofactor $$m = \lceil\frac{t}{p}\rceil$$ to construct $$n = v_p(N!)$$ factors $$mp \ge t$$ that are as small as possible (given that they are divisible by $p$).  To improve efficiency, the algorithm subdivides the interval $$[\sqrt{t},N]$$ into regions within which the values of $m$ and $n$ are constant and simply counts the primes in each region (for small regions it uses `primesieve` to enuemrate primes, for large regions $$[a,b]$$ it uses `primecount` to compute $$\pi(b)-\pi(a-1))$$.

For primes below $$\sqrt{t}$$, the standard greedy approach is modified to only consider cofactors $$m$$ that are $$\sqrt{t}$$-smooth and bounded by $$t^\frac{5}{8}$$ using a precomputed table of factorizations of all such $$m$$ (even for $$10^{12}$$ this takes less than 300MB of memory).  The time spent on this phase of the algorithm is negligible.

To verify the Erdos-Guy-Selfridge conjecture for $$N$$ in $$[10^6,10^{12}]$$, after installing primesieve and primecount and compiling `egs.c` using `build.sh`, you can type

`./egs -f -h hint_1e6_1e12_heuristic_fast.txt`

this should take a minute or two.  If you only want to verify the range $$[10^6,10^{11}]$$ use

`./egs -h -f -h hint_1e6_1e11_heuristic_fast.txt`

which should take 15-20 seconds.  Note that you need to include the `-f` option to use the fast variant of the greedy algorithm, as the standard greedy algorithm will be used by default.

The *hint-files* contain lists of pairs $$(N,t)$$ such that the algorithm is known to produce at least $$N$$ factors of size at least $$t > N/3$$ on these inputs, such that the range [10^6,10^12] is completely covered by the intervals $$[N,N+(t-\lceil \frac{N}{3}\rceil)]$$ (here we are exploiting $$t(N+1) \ge t(N)+1)$$.

To recreate the file `hint_1e6_1e12_heuristic_fast.txt` (or create a new one for a different range) you can type

`./egs -f -c -h hint_1e6_1e12_heuristic_fast.txt 1e6-1e12

This should take 20-25 minutes (recreating `hint_1e6_1e12_heuristic_fast.txt` takes about 5 minutes).

You can compute lower bounds on $$t(N)$$ for particular $$N$$ by typing `egs N` or `egs -f N` (this will use a heuristic search for the best bound that can be easily found)  or use `egs -o N` or `egs -f -o N` to do an exhaustive search for the best bound that the algorithm can prove for the given $$N$$ (this computation will take longer and use all available cores).

```
$ ./egs 1e6
t(1000000) >= 340790 (heuristic greedy) with (t-ceil(N/3)) = 7456 (0.013s)
$ ./egs -f 1e6
t(1000000) >= 337642 (heuristic fast) with (t-ceil(N/3)) = 4308 (0.001s)
$ ./egs -o 1e6
t(1000000) >= 342303 (exhaustive greedy) with (t-ceil(N/3)) = 8969 (0.343s)
$ ./egs -o -f 1e6
t(1000000) >= 338782 (exhaustive fast) with (t-ceil(N/3)) = 5448 (0.065s)

```

To see all command options, run `egs` without specifying any arguments:
```
$ ./egs
Usage: egs [-v level] [-h filename] [-d filename] [-r] [-c] [-o] [-f] N-range [t]
       -v level      integer verbosity level -1 to 4 (optional, default is 0)
       -h filename   hint-file with records N:t (required if range of N is specified)
       -d filename   output-file to dump factorization to (one factor per line, only valid if t is specified)
       -r            verify factorization (set automatically if dump is specified)
       -c            create hint-file rather than reading it (must be specified in combination with -h)
       -o            use the best t for which the algorithm can prove t(N) >= t (optional)
       -f            use fast version of greedy algorithm
       N-range       integer N or range of integers minN-maxN (required, scientific notation supported)
       t             integer t to use for single N (optional, a good t will be determined if unspecified)
```
