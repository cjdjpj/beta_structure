import scipy
import numpy as np

def T2(a, N):
    """Returns the expected pairwise coalescence time in a Beta coalescent"""
    return np.power(1 + 1 / np.exp2(a - 1) / (a - 1), a) * np.power(N, a - 1) / a / scipy.special.beta(2 - a, a)

def n_beta(a, T2):
    """Returns the N necessary to make a Beta coalescent with the given alpha = a have the specified pairwise coalescence time"""
    """if returns 0, means N is unimportant. if returns inf, means T unattainable with given alpha"""
    return ((T2 * a * scipy.special.beta(2 - a, a)) / ((1 + 1 / (2**(a - 1) * (a - 1)))**a))**(1 / (a - 1))

print(n_beta(1.9,25000000))

