# smoothfac

This repository contains an implementation of the `smoothfac` algorithm for
calculating factorizations of $`N!`$ with large factors. It is a hybrid
algorithm choosing factors with large prime divisors greedily and using linear
programming to handle smooth factors having only small prime divisors.

## Problem

Let $`N`$ be a natural number and let $`t`$ be a fixed threshold. The algorithm
finds a subfactorization of $`N!`$ with admissible factors from the set $`J =
\{ t, t+1, \ldots, N \}`$ such that
```math
\prod_{j \in J} j^{x_j} \Big| N!
```
where the $`x_j`$ are non-negative integers counting how often $`j`$ is included in the
product (the notation $`a|b`$ means $`a`$ divides $`b`$).

The objective is to maximize the total number of factors used, *i.e.* the sum
of the $`x_j`$.

A valid subfactorization of $`N!`$ satisfies
```math
\sum_{j \in J} \nu_p(j) \cdot x_j \leq \nu_p(N!)
```
for each prime $`p \in \Pi_N`$ (where $`\Pi_N`$ is the set of primes up to $`N`$ and
$`\nu_p(j)`$ counts the number of times the prime $`p`$ divides the number $`j`$).

The objective together with the constraints for each prime $`p`$ define an integer
linear program. When relaxing the $`x_j`$ to be non-negative reals, the problem
becomes a linear program.

For further discussion, it is useful to switch matrix notation. Let
$`F`$ be the matrix with elements $`(F)_{ij} = \nu_i(j)`$ for $`i \in \Pi_N`$ and
$`j \in J`$, let $`c`$ be the vector with elements $`(c)_i = \nu_i(N!)`$ for
$`i \in \Pi_N`$ and let $`e`$ be the vector of all ones (with conforming dimension). Then, the problem can be written as
```math
\begin{array}{ll}
\displaystyle\max_{x} &  e^\mathsf{T} x \\
\begin{align*}\,\mathrm{s.t.} \\ \phantom{}\end{align*} &
\begin{align*}
    F x &\leq c \\
      x &\geq 0
\end{align*}
\end{array}
```

From any feasible value $`x`$ of the linear program, one can recover a subfactorization
of $`N!`$ as
```math
\prod_{j \in J} j^{\lfloor x_j \rfloor}
```

## Algorithm

Partition the factors $`J = J_S \cup J_R`$ where $`J_S`$ is the set of
$`\sqrt{N}`$-smooth numbers (i.e. all prime divisors of $`j \in J_S`$ are
smaller than or equal to $`\sqrt{N}`$), and $`J_R`$ contains all factors with
a prime divisor larger than $`\sqrt{N}`$.

Partition the prime divisors in small primes $`\Pi_S = \Pi_{\sqrt N}`$ and
large primes $`\Pi_{L} = \Pi_N \backslash \Pi_{\sqrt N}`$.

This allows to rewrite the inequality $`Fx \leq c`$ in block matrix form as
```math
\begin{bmatrix}
  F_{S,S} & F_{S,R} \\
     0    & F_{L,R}
\end{bmatrix}
\cdot
\begin{bmatrix} x_S \\ x_R \end{bmatrix}
\leq
\begin{bmatrix} c_S \\ c_L \end{bmatrix}
```

The first step of the algorithm is to fix values for $`x_R`$ greedily as
follows. Partition the non-smooth factors $`J_R`$ according to their largest
prime divisor
```math
J_R = \bigcup_{p \in \Pi_L} J_p \qquad \mathrm{with} \qquad J_p = \{ j \,|\, \nu_p(j) = 1 \}
```
Choose for each $`J_p`$ its *smallest* element and assign it the full weight
$`c_p`$ of the right hand side
```math
x_j =
\begin{cases}
  \nu_p(N!) & \text{if} \enspace j = p \cdot \lceil t/p \rceil \\
  0         & \text{otherwise,}
\end{cases}
\qquad j \in J_p, \; p \in \Pi_L
```

Fixing $`x_R`$ as above allows to *deflate* the problem and work with a reduced
problem
```math
\begin{array}{ll}
\displaystyle\max_{x_S} &  e^\mathsf{T} x_S \\
\begin{align*}\,\mathrm{s.t.} \\ \phantom{}\end{align*} &
\begin{align*}
    F_{S,S} \cdot x_S &\leq d_S \\
                  x_S &\geq 0
\end{align*}
\end{array}
```
where the deflated right hand side $`d_S = c_S - F_{S,R} \cdot x_R \geq 0`$
[proof needed]. This linear program has only $`\pi(\sqrt{N})`$ constraints
which is a significant reduction from the original $`\pi(N)`$ constraints.

The reduced linear program can be efficiently solved using a sifting strategy.
Start solving the problem with a subset of columns $`W`$ (the *working set*)
containing the $`2\sqrt{N}`$ smallest elements of $`J_S`$. After each solve,
scan for columns not included in the working set having positive reduced cost
```math
  e - w^🞰_W \cdot F_{S,S \backslash W}
```
where $`w^🞰_W`$ are the optimal dual variables from the problem with working
set $`W`$. Scanning is done from smallest to largest element in $`J_S \backslash W`$.
Whenever a batch of 200 improving columns is found, reoptimize the problem
with the columns added to the working set.

Let $`x^🞰_S`$ be an optimal solution of the reduced linear program. The final
step of the algorithm is to greedily use the remaining prime divisors
$`d_S - F_{S,S} \cdot \lfloor x^🞰_S \rfloor`$ to find additional factors
$`j \geq t`$. To this end, pick the largest available prime divisor $`p`$ and
scan if any of the factors $`k p \cdot \lceil t/p \rceil`$, for $`k=1, 2, \ldots`$
divides the product of the remaining prime divisors. If yes, add the factor,
remove the divisors and reiterate; otherwise stop.


## Implementation

The implementation attempts to be space efficient. The bulk of memory is used
for one list `f[]` of length $`N`$. The elements stored in `f[]` are in the
interval $`[-\sqrt{N}, \sqrt{N}]`$ which allows to use 32-bit signed integers for the desired range $`N \leq 10^{11}`$. All other lists and dictionaries have
size $`O(\sqrt{N})`$, except the explicit representation of the linear program
which is (roughly) $`O(\sqrt{N} \cdot \log{N})`$ (assuming the size of the
working set is a bounded multiple of $`\sqrt{N}`$).

The list `f[]` is used to represent the matrices $`F_{S,S}`$, $`F_{S,R}`$ and
to store the counts $`\nu_p(N!)`$ of the large primes $`p \in \Pi_L`$. The
encoding is as follows:
```math
\texttt{f[j]} =
\begin{cases}
  \text{smallest} \, p | j & \text{if } j \in J_S \\
  -\nu_j(N!)               & \text{if } j \in \Pi_L \\
  -\lceil t/p \rceil       & \text{if } j = p \cdot \lceil t/p \rceil \text{, } p \in \Pi_L \text{, } p < t \\
  0                        & \text{otherwise}
\end{cases}
```
The sign of `f[j]` encodes whether $`j`$ is smooth or non-smooth. To determine
whether $`j \in \Pi_L`$ for $`j \geq t`$, a trial division `j/-f[j]` is
required.

The columns of $`F_{S,S}`$ can be recovered by repeated division  
`j ← j / f[j]`.  

The same is possible for the columns of $`F_{S,R}`$ for which the corresponding
variable in $`x_R`$ is non-zero. For a non-smooth composite $`j`$, `f[j]`
contains the ($`\sqrt{N}`$-smooth) value $`-\lceil t/p \rceil`$. Thus, the
repeated division strategy can be started at `-f[j]`.


