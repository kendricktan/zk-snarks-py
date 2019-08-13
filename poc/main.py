"""
https://arxiv.org/abs/1906.07221
"""

import math

from random import randint
from typing import Tuple, List, Union

from py_ecc import bn128
from py_ecc.bn128.bn128_field_elements import inv
from py_ecc.bn128 import FQ2, FQ, add, multiply, double, twist, pairing
from py_ecc.bn128.bn128_pairing import miller_loop, cast_point_to_fq12

# miller_loop        :: Point2D[FQ12] -> Point2D[FQ12] -> FQ12
# twist              :: Point2D[FQ2]  -> Point2D[FQ12]
# cast_point_to_fq12 :: Point2D[FQ]   -> Point2D[FQ12]


# Field order
N: int = bn128.curve_order
P: int = bn128.field_modulus


# Helper functions
def randsn():
    return randint(1, N - 1)


def randsp():
    return randint(1, P - 1)


def mulmodn(x, y):
    return (x * y) % N


def mulmodp(x, y):
    return (x * y) % P


""" 
    Section 3.4 (Page 16)

    Simple restriction for a single coefficient polynomial
"""

alpha = randsp()
s = randsp()
coefficient = randsp()

# Verifier
g_s = multiply(bn128.G1, s)
g_alpha_s = multiply(bn128.G1, alpha * s)

# Prover (doesn't have alpha)
g_s_c = multiply(g_s, coefficient)
g_alpha_s_c = multiply(g_alpha_s, coefficient)

assert multiply(g_s_c, alpha) == g_alpha_s_c


"""
    Section 3.4 (Page 17)

    Restriction for polynomials with multiple coefficients
"""
polynomial_degree = 4

G = 446827
alpha = randsp()
s = randsp()

coefficients = [randint(1, 42) for i in range(polynomial_degree)]

g_s = [
    pow(G, s, P) for i in range(polynomial_degree)
]
g_alpha_s = [
    pow(G, alpha * s, P) for i in range(polynomial_degree)
]

# Evaluate polynomial with provided powers of s
g_p_s_prime = [
    pow(g_s[i], coefficients[i], P) for i in range(polynomial_degree)
]
g_p_s = g_p_s_prime[0]

for c in g_p_s_prime[1:]:
    g_p_s = mulmodp(g_p_s, c)

g_alpha_p_s_prime = [
    pow(g_alpha_s[i], coefficients[i], P) for i in range(polynomial_degree)
]
g_alpha_p_s = g_alpha_p_s_prime[0]

for c in g_alpha_p_s_prime[1:]:
    g_alpha_p_s = mulmodp(g_alpha_p_s, c)

assert pow(g_p_s, alpha, P) == g_alpha_p_s