"""
https://arxiv.org/abs/1906.07221
"""

import math

from functools import reduce
from random import randint
from typing import Tuple, List, Union

from py_ecc import bn128
from py_ecc.bn128.bn128_field_elements import inv
from py_ecc.bn128 import FQ2, FQ, add, multiply, double, twist, pairing
from py_ecc.bn128.bn128_pairing import miller_loop, cast_point_to_fq12


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
    pow(G, pow(s, i, P), P) for i in range(polynomial_degree)
]
g_alpha_s = [
    pow(G, alpha * pow(s, i, P), P) for i in range(polynomial_degree)
]

# Evaluate polynomial with provided powers of s
g_p_s = reduce(
    lambda acc, x: mulmodp(acc, x),
    [pow(g_s[i], coefficients[i], P) for i in range(polynomial_degree)]
)

g_alpha_p_s = reduce(
    lambda acc, x: mulmodp(acc, x),
    [pow(g_alpha_s[i], coefficients[i], P) for i in range(polynomial_degree)]
)

assert pow(g_p_s, alpha, P) == g_alpha_p_s


"""
    Section 3.7 (Page 24)
"""

# Coefficients (a.k.a toxic wastes)
# p(x)/t(x) = h(x)
# (they're in reverse order due to left-to-right as opposed to right-to-left)

polynomial_degree = 4

# p(x) = 1x^3 + 3x^2 + 2x + 0
p_coefficients = [1, 3, 2, 0][::-1]

# t(x) = 0x^3 + 1x^2 - 3x + 2
t_coefficients = [0, 1, 3, 2][::-1]

# h(x) = 0x^3 + 0x^2 - 1x + 0
h_coefficients = [0, 0, 1, 0][::-1]

s = 4 # randsp()
alpha = randsp()

### Setup
g_s_i = [multiply(bn128.G2, pow(s, i, P)) for i in range(polynomial_degree)]
g_alpha_s_i = list(map(lambda x: multiply(x, alpha), g_s_i))

g_alpha = multiply(bn128.G1, alpha)
g_t_s = reduce(
    lambda acc, x: add(acc, x),
    [multiply(bn128.G1, t_coefficients[i] * pow(s, i, P) % P) for i in range(polynomial_degree)]
)

### Proving
g_p_s = reduce(
    lambda acc, x: add(acc, x),
    [multiply(g_s_i[i], p_coefficients[i]) for i in range(polynomial_degree)]
)

g_h_s = reduce(
    lambda acc, x: add(acc, x),
    [multiply(g_s_i[i], h_coefficients[i]) for i in range(polynomial_degree)]
)

g_alpha_p_s = reduce(
    lambda acc, x: add(acc, x),
    [multiply(g_alpha_s_i[i], p_coefficients[i]) for i in range(polynomial_degree)]
)


# Shift polynomial
theta = randsp()

g_theta_p_s = multiply(g_p_s, theta)
g_theta_h_s = multiply(g_h_s, theta)
g_theta_alpha_p_s = multiply(g_alpha_p_s, theta)


### Verification
# add(g^a, g^2a) == multiply(g^a, 3)
#
# miller_loop        :: Point2D[FQ12] -> Point2D[FQ12] -> FQ12
# twist              :: Point2D[FQ2]  -> Point2D[FQ12]
# cast_point_to_fq12 :: Point2D[FQ]   -> Point2D[FQ12]

# a = randsp()

# # g^a
# g_a = multiply(bn128.G1, a)

# # g^a^2 = g^2a
# g_2_a = multiply(g_a, 2)

# # g^2a + g^a = g^3a
# assert add(g_a, g_2_a) == multiply(g_a, 3)

e_gpprime_g = pairing(
    g_theta_alpha_p_s,
    bn128.G1
)

e_gp_galpha = pairing(
    g_theta_p_s,
    g_alpha
)

print(e_gpprime_g)
print(e_gp_galpha)
passed_poly_restrict = e_gpprime_g == e_gp_galpha
print(f'Polynomial restrictions passed: {passed_poly_restrict}')

e_gp_g = pairing(
    g_theta_p_s,
    bn128.G1
)

e_gts_gh = pairing(
    g_theta_h_s,
    g_t_s,
)

passed_poly_cofactors = e_gp_g == e_gts_gh
print(e_gp_g)
print(e_gts_gh)
print(f'Polynomial cofactors passed: {passed_poly_cofactors}')