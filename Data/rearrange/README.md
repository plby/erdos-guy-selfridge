# Data related to the rearrangement section

This directory contains 5 files related to Section "6. Rearranging the
standard factorization":

* `up_to_2.txt` contains values of $t_2(N)$ for $N\le 100000$.  Note
  that the program `src/rearrange/rearrange_2.cc` can compute values
  of $t_2(N)$ very quickly, far beyond this range.

* `up_to_3.txt` contains values of $t_{2,3}(N)$ for $N\le 50000$.
  Note that the largest value for which $t_{2,3}(N) \ge N/4$ is
  $N=26244$.

* `up_to_5.txt` contains values of $t_{2,3,5}(N)$ for $N\le 50000$.
  We plot these numbers in the paper, but do not use the data
  otherwise, because the limiting value of $t_{2,3,5}(N)/N$ appears to
  be below $2/7$.

* `up_to_7.txt` contains values of $t_{2,3,5,7}(N)$ for $N\le 50000$.
  Note that $t_{2,3,5,7}(N) \ge \lfloor 2N/7\rfloor$ in this range.

* `up_to_7.ladder.txt` contains sparse values of $t_{2,3,5,7}(N)$ that
  prove $t_{2,3,5,7}(N) \ge \lfloor 2N/7\rfloor$ for $56<N\le12\times
  10^6$.  The fact that $t_{2,3,5,7}(N+1) \ge t_{2,3,5,7}(N)$ is used
  to minimize the number of terms in the "ladder".
