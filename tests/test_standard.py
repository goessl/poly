from poly import *
import numpy as np
from random import randint, random
from functools import reduce
from itertools import count, islice
from fractions import Fraction
from vector import vecrand
import sympy as sp
from sympy.abc import x as spx
import pytest



#creation
def test_polyzero():
    assert polyzero == ()

def test_polyone():
    assert polyone == (1,)

def test_polyx():
    assert polyx == (0, 1)

def test_polymono():
    assert polymono(3, 5) == (0, 0, 0, 5)
    assert polymono(-1, 25) == polyzero

def test_polymonos():
    for n, xn in enumerate(islice(polymonos(), 10)):
        assert polymono(n) == xn

def test_polyrand():
    for n in range(10):
        p = polyrand(n)
        assert isinstance(p, tuple) and polydeg(p)==n
        assert all(isinstance(ai, float) and 0<=ai<1 for ai in p)

def test_polyrandn():
    for _ in range(10):
        n = randint(0, 10)
        p = polyrandn(n)
        assert isinstance(p, tuple) and polydeg(p)==n
        assert all(isinstance(ai, float) for ai in p)
    assert polyrandn(-1, normed=False) == polyzero

def test_polyfromroots():
    for _ in range(100):
        n = randint(0, 10)
        r = vecrand(n)
        assert np.allclose(polyfromroots(*r),
                           np.polynomial.polynomial.polyfromroots(r))
    assert polyfromroots() == polyone


#utility
def test_polydeg():
    assert polydeg(polyzero) == -1
    assert polydeg(polyone) == 0
    assert polydeg(polyx) == 1
    assert polydeg(polymono(3)) == 3

def test_polyeq():
    assert polyeq((), ())
    assert polyeq((0,), ())
    assert polyeq((0, 0), (0,))
    assert polyeq((1,), (1, 0))
    assert not polyeq((1,), ())
    assert not polyeq((1, 2, 3), (1, 2, 4))

def test_polytrim():
    assert polytrim(()) == ()
    assert polytrim((0,)) == ()
    assert polytrim((1, 0)) == (1,)

def test_polyround():
    assert polyround(()) == ()
    assert polyround((1.1, 2.2)) == (1, 2)
    assert polyround((1.12, 2.23), ndigits=1) == (1.1, 2.2)


#evaluation
def test_polyval():
    #compare with numpy
    for _ in range(1000):
        p, x = polyrand(randint(1, 10)), random()
        actual = np.polynomial.polynomial.polyval(x, p)
        for method in ('naive', 'iterative', 'horner'):
            assert np.isclose(polyval(p, x, method), actual)
    
    #empty polynomial
    assert polyval(polyzero, x, 'naive') \
            == polyval(polyzero, x, 'iterative') \
            == polyval(polyzero, x, 'horner') == 0
    
    #type consistency
    p, x = polyzero, Fraction(5, 2)
    for method in ('naive', 'iterative', 'horner'):
        assert polyval(p, x, method) == 0 \
        and isinstance(polyval(p, x, method), Fraction)
    p, y = (1, 2, 3), Fraction(99, 4)
    for method in ('naive', 'iterative', 'horner'):
        assert polyval(p, x, method) == y \
                and isinstance(polyval(p, x, method), Fraction)
    
    #method error handling
    with pytest.raises(ValueError, match='Invalid method'):
        polyval((1, 2, 3), 4, method='I dont want to do this anymore')

#polyvals gets tested with polyval_iterative

def test_polyvalzero():
    assert polyvalzero(polyzero) == 0
    assert polyvalzero(polyone) == 1
    assert polyvalzero(polyx) == 0
    assert polyvalzero((5, 4, 3)) == 5

def test_polycom():
    def nppolycom(p, q):
        p = np.polynomial.polynomial.Polynomial(p)
        q = np.polynomial.polynomial.Polynomial(q)
        r = p(q)
        return r.coef
    
    for _ in range(1000):
        p, q = polyrand(randint(1, 10)), polyrand(randint(1, 10))
        prediction0 = polycom(p, q, 'naive')
        prediction1 = polycom(p, q, 'iterative')
        prediction2 = polycom(p, q, 'horner')
        actual = nppolycom(p, q)
        assert np.allclose(prediction0, actual)
        assert np.allclose(prediction1, actual)
        assert np.allclose(prediction2, actual)
    
    for method in {'naive', 'iterative', 'horner'}:
        assert polycom(polyzero, polyzero, method) == polyzero
        assert polycom(polyzero, polyone, method) == polyzero
        assert polycom((4, 5, 6), polyzero, method) == (4,)
    
    with pytest.raises(ValueError, match='Invalid method'):
        polycom((1, 2, 3), (4, 5, 6), method='I dont want to do this anymore')

def test_polyshift():
    assert polyshift(polyzero, 5) == polyzero
    assert polyshift(polyone, 5) == polyone
    assert polyshift((1, 2, 3), 5) == polycom((1, 2, 3), polyaddc(polyx, -5))

def test_polyscale():
    assert polyscale(polyzero, 5) == polyzero
    assert polyscale(polyone, 5) == polyone
    assert polyscale((1, 2, 3), 5) == polycom((1, 2, 3), polyscalarmul(5, polyx))


#arithmetic
def test_polypos():
    assert polypos((+1, -2, +3)) == (+1, -2, +3)

def test_polyneg():
    assert polyneg((+1, -2, +3)) == (-1, +2, -3)

def test_polyadd():
    assert polyadd() == ()
    assert polyadd((1, 2)) == (1, 2)
    assert polyadd((1, 2, 3), (4, 5)) == (5, 7, 3)

def test_polyaddc():
    assert polyaddc((), 2, 3) == (0, 0, 0, 2)
    assert polyaddc((1, 2), 4, 5) == (1, 2, 0, 0, 0, 4)
    assert polyaddc((1, 2, 3, 4, 5), 5, 2) == (1, 2, 8, 4, 5)

def test_polysub():
    assert polysub((), ()) == ()
    assert polysub((1,), ()) == (1,)
    assert polysub((), (1,)) == (-1,)
    assert polysub((1, 2, 3), (4, 5)) == (-3, -3, 3)

def test_polyscalarmul():
    assert polyscalarmul(5, polyzero) == polyzero
    assert polyscalarmul(5, (1, 2, 3)) == (5, 10, 15)

def test_polyscalartruediv():
    assert polyscalartruediv((), 1) == ()
    assert polyscalartruediv((4,), 2) == (2,)

def test_polyscalarfloordiv():
    assert polyscalarfloordiv((), 1) == ()
    assert polyscalarfloordiv((3,), 2) == (1,)

def test_polyscalarmod():
    assert polyscalarmod((), 2) == ()
    assert polyscalarmod((3,), 2) == (1,)

def test_polydivmod():
    p, a = (1, 2, 3), 2
    assert polyscalardivmod(p, a) == (polyscalarfloordiv(p, a), polyscalarmod(p, a))

def test_polymul():
    for method in {'naive', 'karatsuba'}:
        for _ in range(1000):
            ps = [polyrand(randint(0, 10)) for _ in range(randint(1, 4))]
            prediction = polymul(*ps, method=method)
            actual = reduce(np.polynomial.polynomial.polymul, ps)
            assert np.allclose(prediction, actual)
        p = polyrand(randint(1, 10))
        assert polymul(method=method) == polyone
        assert polymul(p, method=method) == p
        assert polymul(p, polyzero, method=method) == polyzero
        assert polymul(p, polyone, method=method) == p
    
    with pytest.raises(ValueError, match='Invalid method'):
        polymul((1, 2, 3), (4, 5, 6), method='I dont want to do this anymore')

def test_polymulx():
    for _ in range(1000):
        p = polyrand(randint(1, 10))
        assert np.allclose(polymulx(p), np.polynomial.polynomial.polymulx(p))

def test_polypow():
    for method in {'naive', 'binary'}:
        for _ in range(1000):
            p, n = polyrand(randint(1, 4)), randint(0, 3)
            assert np.allclose(polypow(p, n, method),
                    np.polynomial.polynomial.polypow(p, n))
        assert polyeq(polypow(p, 0, method), polyone)
        assert polyeq(polypow(p, 1, method), p)
    
    with pytest.raises(ValueError, match='Invalid method'):
        polypow((1, 2, 3), 4, method='I dont want to do this anymore')

#def test_polydiv():
#    for _ in range(1000):
#        p, q = polyrand(randint(1, 10)), polyrand(randint(1, 10))
#        prediction = polydiv(tuple(p), tuple(q))
#        actual = np.polynomial.polynomial.polydiv(p, q)
#        assert np.allclose(prediction[0], actual[0]) and np.allclose(prediction[1], actual[1])
#    assert all(map(polyeq, polydiv(polyzero, p), (polyzero, polyzero)))
#    assert all(map(polyeq, polydiv(p, polyone), (p, polyzero)))


def test_polyder():
    for _ in range(1000):
        p = polyrand(randint(0, 10))
        k = randint(1, 10)
        assert np.allclose(polyder(p, k),
                np.polynomial.polynomial.polyder(p, k))
    assert polyder(polyzero) == polyzero
    assert polyder(polyone) == polyzero

def test_polyantider():
    for _ in range(1000):
        p = polyrand(randint(0, 10))
        b, c = random(), random()
        assert np.allclose(polyantider(p, b=b, c=c),
                np.polynomial.polynomial.polyint(p, k=c, lbnd=b))
    assert polytrim(polyantider(polyzero)) == polyzero
    assert polyantider(polyzero, c=1) == polyone
    assert polyantider(polyone) == polyx


#sympy
def test_polysympify():
    assert polysympify((1, 2, 3)) == sp.Poly(1+2*spx+3*spx**2)
    assert polysympify(polyzero) == 0

def test_polyunsympify():
    assert polyunsympify(sp.Poly(1+2*spx+3*spx**2)) == (1, 2, 3)
