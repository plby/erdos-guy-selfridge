# you need to install primesieve version 12.6 or later and primecount version 7.14 or later
gcc -O3 -fopenmp -o egs egs.c -lprimesieve -lprimecount -lstdc++ -lm
