import math

from random import randint
from typing import Tuple, List, Union

from py_ecc import bn128
from py_ecc.bn128.bn128_field_elements import inv
from py_ecc.bn128 import FQ2, FQ, add, multiply, double, twist, pairing
from py_ecc.bn128.bn128_pairing import miller_loop, cast_point_to_fq12


# Types
Point = Tuple[int, int]
PointQ2 = Tuple[Point, Point]  # Quadratic Point with degree 2
PointQ12 = List[int]
AnyPoint = Union[Point, PointQ2, PointQ12]
Scalar = int


# Field order
N: int = bn128.curve_order
P: int = bn128.field_modulus


# Helper functions
def ascoeff(x):
    return x.coeffs if isinstance(x, bn128.FQ2) or isinstance(x, bn128.FQ12) else x


def asint(x):
    return x.n if isinstance(x, bn128.FQ) else x


def fqToPoint(x: AnyPoint):
    if isinstance(x, tuple):
        if isinstance(x[0], bn128.FQ):
            return (asint(x[0]), asint(x[1]))
        if isinstance(x[0], FQ2):
            return (ascoeff(x[0]), ascoeff(x[1]))
        elif isinstance(x[0], bn128.FQ12):
            return ascoeff(x)

    t = type(x)
    raise Exception(f'fqToPoint cant convert {t}')


def randsn():
    return randint(1, N - 1)


def randsp():
    return randint(1, P - 1)


def sbmul(s):
    return bn128.multiply(G, asint(s))


def addmodn(x, y):
    return (x + y) % N


def addmodp(x, y):
    return (x + y) % P


def mulmodn(x, y):
    return (x * y) % N


def mulmodp(x, y):
    return (x * y) % P


def submodn(x, y):
    return (x - y) % N


def submodp(x, y):
    return (x - y) % P


def invmodn(x):
    return inv(x, N)


def invmodp(x):
    return inv(x, P)


def negp(x):
    return (x[0], -x[1])


# Pairing tests
G1 = bn128.G1
G2 = bn128.G2

# Polynomials


def p_x(x):
    return pow(x, 3, N) - (3 * pow(x, 2, N)) + (2 * x)


def t_x(x):
    return (x - 1) * (x - 2)


def h_x(x):
    return p_x(x) / t_x(x)


n = 10
s = randsn()
alpha = randsn()

proving_key = (
    [multiply(G1, pow(s, i, N)) for i in range(n)],
    [multiply(G1, alpha * pow(s, i, N) % N) for i in range(n)]
)

verification_key = (
    multiply(G1, t_x(s)),
    multiply(G1, alpha)
)

# Prover evaluates
# prover_Gp = reduce

# miller_loop        :: Point2D[FQ12] -> Point2D[FQ12] -> FQ12
# twist              :: Point2D[FQ2] -> Point2D[FQ12]
# cast_point_to_fq12 :: Point2D[FQ] -> Point2D[FQ12]
sk = randsn()
pk1 = multiply(G1, sk)
pk2 = multiply(G2, sk)

k1 = miller_loop(
    twist(G2),
    cast_point_to_fq12(pk1)
)

k2 = miller_loop(
    twist(pk2),
    cast_point_to_fq12(G1)
)

verified = k1 == k2
print(f'Verified: {verified}')

# miller_loop(
#     cast_point_to_fq12(proving_key[0][0]),
#     cast_point_to_fq12(proving_key[0][1]),
# )
