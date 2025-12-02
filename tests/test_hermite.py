from poly import *
import numpy as np
from vector import *
from random import gauss, randint
from functools import reduce
from itertools import islice
from fractions import Fraction
from math import sqrt, pi, factorial, isclose
import sympy as sp
from sympy.abc import x as spx



#creation
def test_H0():
    assert H0 == (1,)

def test_H1():
    assert H1 == (0, 2)

def test_H2():
    assert H2 == (-2, 0, 4)

def test_herm():
    for method in {'recursive', 'iterative', 'explicit'}:
        for n in range(25): #numpy fails at n=29
            prediction = herm(n, method)
            actual = np.polynomial.hermite.herm2poly(vecbasis(n))
            assert np.array_equal(prediction, actual)
    assert tuple(islice(herms(), 10)) == tuple(herm(i) for i in range(10))

def test_hermzero():
    assert hermzero == ()

def test_hermone():
    assert hermone == (1,)

def test_hermx():
    assert hermx == (0, Fraction(1, 2))

def test_hermmono():
    for method in {'recursive', 'iterative', 'explicit'}:
        for n in range(10):
            prediction = hermmono(n, method)
            actual = np.polynomial.hermite.poly2herm(polymono(n))
            assert np.allclose(np.array(prediction, dtype=float), actual)
    assert tuple(islice(hermmonos(), 10)) == tuple(hermmono(i) for i in range(10))

def test_hermrand():
    for n in range(-1, 10):
        h = hermrand(n)
        assert isinstance(h, tuple) and hermdeg(h)==n

def test_hermrandn():
    assert hermrandn(-1, normed=False) == ()
    for n in range(10):
        h = hermrandn(n)
        assert isinstance(h, tuple) and hermdeg(h)==n and isclose(hermabs(h), 1)

def test_hermfromroots():
    for _ in range(100):
        r = tuple(randint(-100, +100) for _ in range(randint(1, 10)))
        h = hermfromroots(*r)
        assert all(hermval(h, ri)==0 for ri in r)
    assert hermfromroots() == hermone

#Hilbert space
def test_hermweighti():
    for n in range(20):
        assert hermweighti(n) == 2**n*factorial(n)

def test_hermweightis():
    for n, x in enumerate(islice(hermweightis(), 20)):
        assert hermweighti(n) == x

def test_hermweight():
    for n in range(20):
        assert isclose(hermweight(n), sqrt(pi)*2**n*factorial(n))

def test_hermweights():
    for n, x in enumerate(islice(hermweights(), 20)):
        assert isclose(hermweight(n), x)

def test_hermdoti():
    for _ in range(20):
        a = [randint(-100, +100) for _ in range(randint(0, 10))]
        b = [randint(-100, +100) for _ in range(randint(0, 10))]
        
        assert hermdoti(a, b) == sum(ai*bi * 2**i*factorial(i) for i, (ai, bi) in enumerate(zip(a, b)))


#evaluation
def test_hermval():
    for method in {'naive', 'iterative', 'clenshaw'}:
        for _ in range(100):
            h = np.random.rand(np.random.randint(0, 10)).tolist()
            x = 2*np.random.rand()-1
            
            prediction = hermval(h, x, method)
            actual = np.polynomial.hermite.hermval(x, h if h else [0])
            assert np.isclose(prediction, actual)
        assert hermval(hermzero, x) == 0

def test_hermvalzero():
    assert tuple(islice(hermvalzeros(), 20)) == tuple(hermval(vecbasis(i), 0) for i in range(20))
    for _ in range(100):
        h = np.random.randint(-100, +100, size=np.random.randint(0, 10)).tolist()
        prediction = hermvalzero(h)
        actual = hermval(h, 0)
        assert prediction == actual

#arithmetic
def test_hermpos():
    assert hermpos((+1, -2, +3)) == (+1, -2, +3)

def test_hermneg():
    assert hermneg((+1, -2, +3)) == (-1, +2, -3)

def test_hermaddc():
    assert hermaddc((), 2, 3) == (0, 0, 0, 2)
    assert hermaddc((1, 2), 4, 5) == (1, 2, 0, 0, 0, 4)
    assert hermaddc((1, 2, 3, 4, 5), 5, 2) == (1, 2, 8, 4, 5)

def test_hermadd():
    assert hermadd() == ()
    assert hermadd((1, 2)) == (1, 2)
    assert hermadd((1, 2, 3), (4, 5)) == (5, 7, 3)

def test_hermsub():
    assert hermsub((), ()) == ()
    assert hermsub((1,), ()) == (1,)
    assert hermsub((), (1,)) == (-1,)
    assert hermsub((1, 2, 3), (4, 5)) == (-3, -3, 3)

def test_hermscalarmul():
    assert hermscalarmul(5, hermzero) == hermzero
    assert hermscalarmul(5, (1, 2, 3)) == (5, 10, 15)

def test_hermscalartruediv():
    assert hermscalartruediv((), 1) == ()
    assert hermscalartruediv((4,), 2) == (2,)

def test_hermscalarfloordiv():
    assert hermscalarfloordiv((), 1) == ()
    assert hermscalarfloordiv((3,), 2) == (1,)

def test_hermscalarmod():
    assert hermscalarmod((), 2) == ()
    assert hermscalarmod((3,), 2) == (1,)

def test_hermdivmod():
    h, a = (1, 2, 3), 2
    hermscalardivmod(h, a) == (hermscalartruediv(h, a), hermscalarmod(h, a))

def test_hermmul():
    assert hermmul() == hermone
    assert hermmul(hermzero) == hermzero
    assert hermmul(hermone) == hermone
    assert hermmul(hermx) == hermx
    
    for method in {'naive'}:
        for _ in range(100):
            g = np.random.rand(np.random.randint(0, 5))
            h = np.random.rand(np.random.randint(0, 5))
            
            prediction = hermmul(g, h, method=method)
            actual = np.polynomial.hermite.hermmul(g, h) if g.size and h.size else []
            assert np.allclose(prediction, actual)

def test_hermmulx():
    assert hermmulx(hermzero) == (0,)
    for _ in range(100):
        h = np.random.rand(np.random.randint(1, 5))
        
        prediction = hermmulx(h)
        actual = np.polynomial.hermite.hermmulx(h)
        assert np.allclose(prediction, actual)

def test_hermmulHn():
    assert hermmulHn(hermzero, 5) == hermzero
    for _ in range(100):
        h = np.random.rand(np.random.randint(1, 5))
        n = randint(0, 20)
        
        assert np.allclose(hermmulHn(h, n), hermmul(h, vecbasis(n)))

def test_hermpow():
    for method in {'naive', 'binary'}:
        assert hermpow(hermzero, 0) == hermone
        assert hermpow(hermzero, 1) == hermzero
        assert hermpow(hermzero, 2) == hermzero
        assert hermpow(hermone, 0) == hermone
        assert hermpow(hermone, 1) == hermone
        assert hermpow(hermone, 2) == hermone
        
        for _ in range(100):
            h = np.random.rand(np.random.randint(0, 5))
            n = np.random.randint(0, 5)
            
            prediction = hermpow(tuple(h), n, method)
            actual = np.polynomial.hermite.hermpow(h, n) if h.size else []
            assert np.allclose(prediction, actual)

def test_hermmulpow():
    pass


#calculus
def test_hermder():
    for _ in range(100):
        n = randint(0, 10)
        k = randint(0, 10)
        h = np.random.rand(n)
        prediction = hermder(h, k)
        actual = np.polynomial.hermite.hermder(h if len(h) else [], k)
        assert np.allclose((prediction if len(prediction) else [0]), actual)

def test_hermantider():
    for _ in range(100):
        h = np.random.rand(np.random.randint(0, 10))
        
        prediction = hermantider(h)
        actual = np.polynomial.hermite.hermint(h) if h.size else [0]
        actual = hermaddc(actual, hermvalzero(prediction))
        assert np.allclose(prediction, actual)

#sympy
def test_hermsympify():
    assert hermsympify((1, 2, 3)) == sp.Poly.from_list([12, 4, -5], spx)




















def test_polymono():
    assert polymono(3, 5) == (0, 0, 0, 5)

def test_herm2poly():
    for method in {'naive', 'iterative', 'clenshaw'}:
        assert herm2poly(hermzero, method) == polyzero
        for _ in range(100):
            h = np.random.rand(np.random.randint(0, 10)).tolist()
            
            prediction = herm2poly(h, method)
            actual = np.polynomial.hermite.herm2poly(h) if h else []
            assert np.allclose(prediction, actual)

"""
def test_hermmono():
    for method in {'recursive', 'iterative', 'explicit'}:
        for n in range(25):
            prediction = hermmono(n, method)
            assert herm2poly(prediction) == vecbasis(n)

def test_poly2herm():
    for method in {'naive', 'horner'}:
        assert poly2herm((), method) == ()
        for _ in range(100):
            p = np.random.rand(np.random.randint(0, 10)).tolist()
            
            prediction = poly2herm(p, method)
            actual = np.polynomial.hermite.poly2herm(p) if p else []
            assert np.allclose(prediction, actual)









def test_hermintegral():
    for _ in range(100):
        h = np.random.randint(-10, +10, np.random.randint(0, 5))
        a = np.random.randint(1, 10)
        prediction = hermintegral(h, a) * sp.sqrt(sp.pi/a)
        actual = sp.integrate(sp.exp(-a*spx**2)*hermsympify(h), (spx, -sp.oo, +sp.oo)).expand()
        assert prediction == actual

def test_hermpolyintegral():
    for _ in range(100):
        p = np.random.randint(-10, +10, np.random.randint(0, 10))
        a = np.random.randint(1, 10)
        prediction = hermpolyintegral(p, a)
        actual = hermintegral(poly2herm(p), a)
        assert prediction == actual
"""
