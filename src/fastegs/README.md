This directory contains a fast, memory efficient implemention of (a variant of) the greedy algorithm for constructing factorizations of $$N!$$ consisting of integers greater than $$t \ge \frac{N}{3}$$, with the goal of obtaining at least $$N$$ factors (thereby proving $$t(N) \ge t$$).  It is specifically optimized for large values of $$N$$, say $$N \ge 10^6$$; for smaller values of $$N$$ the bounds on $$t$$ it produces will typically not be as good as those given by the standard greedy algorithm, but it is much faster when $$N$$ is large, and it uses only $$\tilde O(N^{\frac{5}{8}})$$ bits of memory (the 5/8 can be replaced by any exponent strictly above 1/2).

This implemention depends on Kim Walisch's [primesieve](https://github.com/kimwalisch/primesieve) and [primecount](https://github.com/kimwalisch/primecount) libraries for enumerating and counting primes, which must be installed before you compile the code in this directory.

These two libraries are used to efficiently handle factors divisible by primes in the range $$[\sqrt{t},N]$$. In this range the greedy algorithm can always use the optimal cofactor $$m = \lceil\frac{t}{p}\rceil$$ to construct $$n = v_p(N!)$$ factors $$mp \ge t$$ that are as small as possible (given that they are divisible by $p$).  To improve efficiency, the algorithm subdivides the interval $$[\sqrt{t},N]$$ into regions within which the values of $m$ and $n$ are constant and simply counts the primes in each region (for small regions it uses `primesieve` to enuemrate primes, for large regions $$[a,b]$$ it uses `primecount` to compute $$\pi(b)-\pi(a-1))$$.

For primes below $$\sqrt{t}$$, the standard greedy approach is modified to only consider cofactors $$m$$ that are $$\sqrt{t}$$-smooth and bounded by $$t^\frac{5}{8}$$ using a precomputed table of factorizations of all such $$m$$ (even for $$10^{12}$$ this takes less than 300MB of memory).  The time spent on this phase of the algorithm is negligible.

To verify the Erdos-Guy-Selfridge conjecture for $$N$$ in $$[10^6,10^{12}]$$, after installing primesieve and primecount and compiling `egs.c` using `build.sh`, you can type

`./egs -h hint12.txt 1e6-1e12`

this should take a minute or two and use substantially less than 1GB memory.  If you only want to verify the range $$[10^6,10^{11}]$$ use

`./egs -h hint11.txt 1e6-1e11`

which should take 15-30 seconds.

The files `hint11.txt` and `hint12.txt` are *hint-files* containing of a list of pairs $$(N,t)$$ such that the algorithm is known to produce at least $$N$$ factors of size at least $$t > N/3$$ on these inputs, such that the range [10^6,10^12] is completely covered by the intervals $$[N,N+(t-\lceil \frac{N}{3}\rceil)]$$ (here we are exploiting $$t(N+1) \ge t(N)+1)$$.

To create recreate the file `hint12.txt` (or create a new one for a different range) you can type

`./egs -c -h hint12.txt 1e6-1e12

This will take an hour or two (recreating `hint11.txt` only takes 5-10 minutes).

You can compute lower bounds on $$t(N)$$ for particular $$N$$ by typing `egs N` or `egs -o N`.

```
$ ./egs 1e6
t(1000000) >= 337642 with (t-ceil(N/3)) = 4308 (0.001s)
$ ./egs -o 1e6
t(1000000) >= 338782 (optimal for algorithm) with (t-ceil(N/3)) = 5448 (0.080s)

```

The former computes a pretty good lower bound using the same method that is used when computing hint files; adding the `-` option asks for the program to determine the last lower bound that can be determined using this variant of the greedy algorithm, which will typically be slightly better, but take much longer to compute.

To see all command options simply don't specify any arguments:
```
$ ./egs
Usage: egs [-v level] [-h filename] [-c] [-o] N-range [t]
       -v level      integer verbosity level -1 to 3 (optional, default is 0)
       -h filename   hint-file with records N:t (required if range of N is specified)
       -c            create hint-file rather than reading it (must be specified in combination with -h)
       -o            use the best t for which the algorithm can prove t(N) >= t (optional)
       N-range       integer N or range of integers M-N (required, scientific notation supported)
       t             integer t to use for single N (optional, a good t will be determined if unspecified)

```
