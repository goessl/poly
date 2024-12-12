from poly import *
import numpy as np
from functools import reduce



def test_polynpzero():
    assert np.array_equal(polynpzero(1), np.array([[0]]))
    assert np.array_equal(polynpzero(2), np.array([[0], [0]]))
    assert np.array_equal(polynpzero(3), np.array([[0], [0], [0]]))

def test_polynpone():
    assert np.array_equal(polynpone(1), np.array([[1]]))
    assert np.array_equal(polynpone(2), np.array([[1], [1]]))
    assert np.array_equal(polynpone(3), np.array([[1], [1], [1]]))

def test_polynpmono():
    assert np.array_equal(polynpmono(0, d=1), np.array([[1]]))
    assert np.array_equal(polynpmono(0, d=2), np.array([[1], [1]]))
    assert np.array_equal(polynpmono(1, d=1), np.array([[0, 1]]))
    assert np.array_equal(polynpmono(1, d=2), np.array([[0, 1], [0, 1]]))


def test_polynpfromroots():
    #assert np.array_equal(polynpfromroots(()), [1])
    
    for _ in range(100): #1D
        r = np.random.rand(np.random.randint(1, 5))
        prediction = polynpfromroots(r)
        actual = np.polynomial.polynomial.polyfromroots(r)
        assert np.allclose(prediction, actual)
    
    for _ in range(100): #2D
        r = np.random.rand(20, np.random.randint(1, 5))
        prediction = polynpfromroots(r)
        actual = [np.polynomial.polynomial.polyfromroots(ri) for ri in r]
        assert np.allclose(prediction, actual)

'''
def test_polynpval():
    for _ in range(100): #1D(0D)
        p = np.random.rand(np.random.randint(1, 10))
        x = np.random.rand()
        prediction = polynpval(p, x)
        actual = np.polynomial.polynomial.polyval(x, p)
        assert np.isclose(prediction, actual)
    
    for _ in range(100): #1D(1D)
        p = np.random.rand(np.random.randint(1, 10))
        x = np.random.rand(100)
        prediction = polynpval(p, x)
        actual = np.polynomial.polynomial.polyval(x, p)
        assert np.allclose(prediction, actual)
    
    for _ in range(100): #2D(0D)
        p = np.random.rand(20, np.random.randint(1, 10))
        x = np.random.rand()
        prediction = polynpval(p, x)
        actual = [np.polynomial.polynomial.polyval(x, pi) for pi in p]
        assert np.allclose(prediction, actual)


def test_polynpcom():
    for _ in range(100):
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        q = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polynpcom(p, q)
        actual = [np.poly1d(pi[::-1])(np.poly1d(qi[::-1])).c[::-1] for pi, qi in zip(p, q)]
        assert np.allclose(prediction, actual)

def test_polynpmul():
    for _ in range(100):
        n = np.random.randint(1, 4)
        ps = [np.random.normal(size=(20, np.random.randint(1, 10))) for _ in range(n)]
        prediction = polynpmul(*ps)
        actual = [reduce(np.polynomial.polynomial.polymul, rows) for rows in zip(*ps)]
        assert np.allclose(prediction, actual)

def test_polynppow():
    for _ in range(100):
        n = np.random.randint(0, 5)
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polynppow(p, n)
        actual = np.apply_along_axis(
                lambda pi: np.polynomial.polynomial.polypow(pi, n), 1, p)
        assert np.allclose(prediction, actual)


def test_polynpder():
    for _ in range(100): #1D
        p = np.random.rand(np.random.randint(1, 10))
        n = np.random.randint(0, 5)
        prediction = polynpder(p, n)
        actual = np.polynomial.polynomial.polyder(p, n)
        assert np.allclose(prediction, actual)
    
    for _ in range(100): #2D
        p = np.random.rand(20, np.random.randint(1, 10))
        n = np.random.randint(0, 5)
        prediction = polynpder(p, n)
        actual = np.apply_along_axis(
                lambda p: np.polynomial.polynomial.polyder(p, n), 1, p)
        assert np.allclose(prediction, actual)

def test_polynpmder():
    for deg in range(10):
        p = polyrand(deg)
        
        prediction = polynpmder(deg) @ p
        actual = polyder(p)
        assert np.allclose(prediction, actual)

def test_polynpint():
    for _ in range(100): #1D
        p = np.random.rand(np.random.randint(1, 10))
        n = np.random.randint(0, 5)
        c = np.random.rand(n)
        prediction = polynpint(p, n, c)
        actual = np.polynomial.polynomial.polyint(p, n, c)
        assert np.allclose(prediction, actual)
    
    for _ in range(100): #2D
        p = np.random.rand(20, np.random.randint(1, 10))
        n = np.random.randint(0, 5)
        c = np.random.rand(n)
        prediction = polynpint(p, n, c)
        actual = np.apply_along_axis(
                lambda p: np.polynomial.polynomial.polyint(p, n, c), 1, p)
        assert np.allclose(prediction, actual)

def test_polynpmint():
    for deg in range(10):
        p = polyrand(deg)
        
        prediction = polynpmint(deg) @ p
        actual = polyint(p)
        assert np.allclose(prediction, actual)
'''
