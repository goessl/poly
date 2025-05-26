from poly import *
import numpy as np
from random import gauss, randint
from vector import vecbasis, vecrand, veceq, vectrim
from functools import reduce
from fractions import Fraction



def test_polyfromroots():
    for _ in range(1000):
        r = vecrand(randint(1, 10))
        p = polyfromroots(*r)
        assert all(np.isclose(polyval(p, ri), 0) for ri in r)
        assert not np.allclose(p, 0)


def test_polyval():
    #compare with numpy
    for _ in range(1000):
        p, x = vecrand(randint(1, 10)), gauss()
        assert np.isclose(polyval(p, x, 'naive'),
                np.polynomial.polynomial.polyval(x, p))
        assert np.isclose(polyval(p, x, 'iterative'),
                np.polynomial.polynomial.polyval(x, p))
        assert np.isclose(polyval(p, x, 'horner'),
                np.polynomial.polynomial.polyval(x, p))
    #empty polynomial
    assert polyval(polyzero, x, 'naive') \
            == polyval(polyzero, x, 'iterative') \
            == polyval(polyzero, x, 'horner') == float(0)
    #type consistency
    p, x = polyzero, Fraction(5, 2)
    for method in ('naive', 'iterative', 'horner'):
        assert polyval(p, x, method) == 0 \
        and isinstance(polyval(p, x, method), Fraction)
    p = (1, 2, 3)
    for method in ('naive', 'iterative', 'horner'):
        assert polyval(p, x, method) == Fraction(99, 4) \
                and isinstance(polyval(p, x, method), Fraction)

def test_polycom():
    def nppolycom(p, q):
        p = np.polynomial.polynomial.Polynomial(p)
        q = np.polynomial.polynomial.Polynomial(q)
        r = p(q)
        return r.coef
    
    for _ in range(1000):
        p, q = vecrand(randint(1, 10)), vecrand(randint(1, 10))
        prediction0 = polycom(p, q, 'naive')
        prediction1 = polycom(p, q, 'iterative')
        prediction2 = polycom(p, q, 'horner')
        actual = nppolycom(p, q)
        assert np.allclose(prediction0, actual)
        assert np.allclose(prediction1, actual)
        assert np.allclose(prediction2, actual)
    
    for method in {'naive', 'iterative'}:
        assert polycom(polyzero, polyzero, method) == polyzero
        assert polycom(polyzero, polyone, method) == polyzero
        assert polycom((4, 5, 6), polyzero, method) == (4,)


def test_polymul():
    for method in {'naive', 'karatsuba'}:
        for _ in range(1000):
            ps = [vecrand(randint(1, 10)) for _ in range(randint(1, 4))]
            prediction = polymul(*ps, method=method)
            actual = reduce(np.polynomial.polynomial.polymul, ps)
            assert np.allclose(prediction, actual)
        p = vecrand(randint(1, 10))
        assert polymul(method=method) == polyone
        assert polymul(p, method=method) == p
        assert vectrim(polymul(p, polyzero, method=method)) == polyzero
        assert polymul(p, polyone, method=method) == p

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
