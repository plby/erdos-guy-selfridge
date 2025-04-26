This directory contains a fast, memory efficient implemention of (a variant of) the greedy algorithm for constructing factorizations of $$N!$$ consisting of integers greater than $$t \ge \frac{N}{3}$$, with the goal of obtaining at least $$N$$ factors (thereby proving $$t(N) \ge t$$).  It is specifically optimized for large values of $$N$$, say $$N \ge 10^6$$; for smaller values of $$N$$ the bounds on $$t$$ it produces will typically not be as good as those given by the standard greedy algorithm, but it is much faster when $$N$$ is large, and it uses only $$\tilde O(N^{\frac{5}{8}})$$ bits of memory (the 5/8 can be replaced by any exponent strictly above 1/2).

This implemention depends on Kim Walisch's [primesieve](https://github.com/kimwalisch/primesieve) and [primecount](https://github.com/kimwalisch/primecount) libraries for enumerating and counting primes, which must be installed before you compile the code in this directory.

These two libraries are used to efficiently handle factors divisible by primes in the range $$[\sqrt{t},N]$$. In this range the greedy algorithm can always use the optimal cofactor $$m = \lceil\frac{t}{p}\rceil$$ to construct $$n = v_p(N!)$$ factors $$mp \ge t$$ that are as small as possible (given that they are divisible by $p$).  To improve efficiency, the algorithm subdivides the interval $$[\sqrt{t},N]$$ into regions within which the values of $m$ and $n$ are constant and simply counts the primes in each region (for small regions it uses `primesieve` to enuemrate primes, for large regions $$[a,b]$$ it uses `primecount` to compute $$\pi(b)-\pi(a-1))$$.

For primes below $$\sqrt{t}$$, the standard greedy approach is modified to only consider cofactors $$m$$ that are $$\sqrt{t}$$-smooth and bounded by $$t^\frac{5}{8}$$ using a precomputed table of factorizations of all such $$m$$ (even for $$10^{12}$$ this takes less than 300MB of memory).  The time spent on this phase of the algorithm is negligible.

To verify the Erdos-Guy-Selfridge conjecture for $$N$$ in $$[10^6,10^{12}]$$, after installing primesieve and primecount and compiling `egs.c` using `build.sh`, you can type

`./egs 1000000-1000000000000 hint12.txt`

this should take a minute or two and use substantially less than 1GB memory.  If you only want to verify the range $$[10^6,10^{11}]$$ you can type

`./egs 1000000-100000000000 hint12.txt`

which should take 15-30 seconds.  To see more detail about what is happening you can set the verbosity parameter, e.g.

`./egs 1000000-100000000000 hint12.txt 1`

The file `hint12.txt` consists of a list of pairs (N,t) such that the algorithm is known to produce at least N factors of size at least t > N/3 on these inputs, such that the range [10^6,10^12] is completely covered by the intervals $$[N,N+(t-\lceil \frac{N}{3}\rceil)]$$ (here we are exploiting $$t(N+1) \ge t(N)+1)$$.

To reconstruct the hint-file from scratch, simply run the program without specifying a hint-file and redirect the output, e.g.

`./egs 1000000-100000000000 >hint11.txt`

This takes substantially longer (5-10 minutes for $$[10^6,10^{11}]$$, an hour or so for $$[10^6,10^{12}]$$).
