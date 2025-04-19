from calculations import evaluate
import numpy as np

MIN_A = 180
MAX_A = 190
JUMP_A = 0.2

MIN_K = 290
MAX_K = 300
JUMP_K = 1

# L seemed to always go towards 4.5
L = 4.5

MIN_N_EXP = 10.0
MAX_N_EXP = 11.1
JUMP_N_EXP = 0.001

best_A = None
best_K = None
best_L = None
best_N_exp = MAX_N_EXP

for A in np.arange(MIN_A, MAX_A, JUMP_A):
    for K in np.arange(MIN_K, MAX_K, JUMP_K):
        for N_EXP in np.arange(MIN_N_EXP, best_N_exp, JUMP_N_EXP)[::-1]:
            N = 10**N_EXP
            t = N/3
            try:
                evaluate(t, N, A, K, L)

                if N_EXP < best_N_exp:
                    best_A = A
                    best_K = K
                    best_L = L
                    best_N_exp = N_EXP
                    print(f"New best: N_EXP={N_EXP}, A={A}, K={K}, L={L}")
            except Exception as e:
                break