from poly import *
import numpy as np
from vector import *
from random import gauss, randint
from functools import reduce
from fractions import Fraction



def test_herm():
    for n in range(25): #fails at n=29
        prediction0 = herm(n, 'recursive')
        prediction1 = herm(n, 'iterative')
        prediction2 = herm(n, 'explicit')
        actual = np.polynomial.hermite.herm2poly(vecbasis(n))
        assert np.array_equal(prediction0, actual)
        assert np.array_equal(prediction1, actual)
        assert np.array_equal(prediction2, actual)

def test_herm2poly():
    for n in range(100):
        h = np.random.rand(np.random.randint(0, 10)).tolist()
        
        prediction0 = herm2poly(h, 'naive')
        prediction1 = herm2poly(h, 'clenshaw')
        actual = np.polynomial.hermite.herm2poly(h) if h else []
        assert np.allclose(prediction0, actual)
        assert np.allclose(prediction1, actual)
    
    assert herm2poly((), 'naive') == herm2poly((), 'clenshaw') == ()

def test_hermmono():
    for n in range(25):
        prediction0 = hermmono(n, 'recursive')
        prediction1 = hermmono(n, 'iterative')
        prediction2 = hermmono(n, 'explicit')
        assert herm2poly(prediction0) == vecbasis(n)
        assert herm2poly(prediction1) == vecbasis(n)
        assert herm2poly(prediction2) == vecbasis(n)

def test_poly2herm():
    pass


def test_hermval():
    for n in range(100):
        h = np.random.rand(np.random.randint(0, 10)).tolist()
        x = 2*np.random.rand()-1
        
        prediction1 = hermval(h, x, 'naive')
        prediction2 = hermval(h, x, 'clenshaw')
        actual = np.polynomial.hermite.hermval(x, h if h else [0])
        assert np.isclose(prediction1, actual)
        assert np.isclose(prediction2, actual)



def test_hermmulx():
    for _ in range(100):
        h = tuple(np.random.rand(np.random.randint(1, 5)))
        prediction = hermmulx(h)
        actual = np.polynomial.hermite.hermmulx(h)
        assert np.allclose(prediction, actual)
    assert hermmulx(()) == (0,)


def test_hermder():
    from numpy.polynomial.hermite import poly2herm, herm2poly
    from numpy.polynomial.polynomial import polyder
    
    for _ in range(100):
        h = np.random.rand(np.random.randint(0, 10)).tolist()
        
        prediction = hermder(h)
        actual = np.polynomial.hermite.hermder(h) if h else []
        assert np.allclose(prediction, actual)
