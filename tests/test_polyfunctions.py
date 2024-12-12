from poly import *
import numpy as np
from random import gauss, randint
from vector import vecbasis, vecrand, veceq, vectrim
from functools import reduce
from fractions import Fraction



def test_polyfromroots():
    for _ in range(1000):
        r = vecrand(randint(1, 10))
        p = polyfromroots(r)
        assert all(np.isclose(polyval(p, ri), 0) for ri in r)
        assert not np.allclose(p, 0)


def test_polyval():
    for _ in range(1000):
        p, x = vecrand(randint(1, 100)), gauss()
        assert np.isclose(polyval(p, x, 'naive'),
                np.polynomial.polynomial.polyval(x, p))
        assert np.isclose(polyval(p, x, 'iterative'),
                np.polynomial.polynomial.polyval(x, p))
        assert np.isclose(polyval(p, x, 'horner'),
                np.polynomial.polynomial.polyval(x, p))
    assert polyval(polyzero, Fraction(5)) == 0 \
            and isinstance(polyval(polyzero, Fraction(5)), Fraction)

def test_polycom():
    for _ in range(1000):
        p, q = vecrand(randint(1, 10)), vecrand(randint(1, 10))
        assert np.allclose(polycom(p, q),
                    np.poly1d(p[::-1])(np.poly1d(q[::-1])).c[::-1])
    assert polycom(polyzero, polyzero) == polyzero
    assert polycom(polyzero, polyone) == polyzero
    assert polycom(polyone, polyzero) == polyone
    assert polycom(polyzero, p) == polyzero
    assert polycom(p, polyzero) == vecbasis(0, p[0])


def test_polymul():
    for _ in range(1000):
        ps = [vecrand(randint(1, 10)) for _ in range(randint(1, 4))]
        assert np.allclose(polymul(*ps),
                reduce(np.polynomial.polynomial.polymul, ps))
    p = vecrand(randint(1, 10))
    assert polymul() == polyone
    assert polymul(p) == p
    assert vectrim(polymul(p, polyzero)) == polyzero
    assert polymul(p, polyone) == p

def test_polymulx():
    for _ in range(1000):
        p = vecrand(randint(1, 10))
        assert np.allclose(polymulx(p), np.polynomial.polynomial.polymulx(p))

def test_polypow():
    for _ in range(1000):
        p, n = vecrand(randint(1, 4)), randint(0, 3)
        assert np.allclose(polypow(p, n),
                np.polynomial.polynomial.polypow(p, n))
    assert veceq(polypow(p, 0), polyone)
    assert veceq(polypow(p, 1), p)

def test_polydiv():
    for _ in range(1000):
        p, q = vecrand(randint(1, 10)), vecrand(randint(1, 10))
        prediction = polydiv(tuple(p), tuple(q))
        actual = np.polynomial.polynomial.polydiv(p, q)
        assert np.allclose(prediction[0], actual[0]) and np.allclose(prediction[1], actual[1])
    assert all(map(veceq, polydiv(polyzero, p), (polyzero, polyzero)))
    assert all(map(veceq, polydiv(p, polyone), (p, polyzero)))


def test_polyder():
    for _ in range(1000):
        p = vecrand(randint(1, 10))
        assert np.allclose(polyder(p),
                np.polynomial.polynomial.polyder(p))
    assert polyder(polyzero) == polyzero
    assert polyder(polyone) == polyzero

def test_polyint():
    for _ in range(1000):
        p, c = vecrand(randint(1, 10)), gauss()
        assert np.allclose(polyint(p, c=c),
                np.polynomial.polynomial.polyint(p, k=c))
    assert vectrim(polyint(polyzero)) == polyzero
    assert polyint(polyzero, c=1) == polyone
    assert polyint(polyone) == polyx
    assert polyint(polyzero, n=2, c=(1, 2)) == (2, 1)
