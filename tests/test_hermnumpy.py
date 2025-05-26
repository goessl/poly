from poly import *
import numpy as np



def test_hermnp2poly():
    for _ in range(100): #1D
        h = np.random.rand(np.random.randint(0, 10))
        
        prediction = hermnp2poly(h)
        actual = herm2poly(h)
        assert np.allclose(prediction, actual)
    
    for _ in range(100): #2D
        hs = np.random.rand(2, np.random.randint(0, 10))
        
        prediction = hermnp2poly(hs)
        actual = [herm2poly(h) for h in hs]
        assert np.allclose(prediction, actual)

def test_hermnp2polyM():
    for _ in range(100):
        deg = np.random.randint(0, 10)
        h = np.random.rand(deg+1)
        
        prediction = hermnp2polyM(deg) @ h
        actual = hermnp2poly(h)
        assert np.allclose(prediction, actual)

def test_polynp2herm():
    for deg in range(20):
        M = hermnp2polyM(deg)
        Minv = polynp2hermM(deg)
        assert np.array_equal(M@Minv, np.eye(deg+1))


def test_hermnpder():
    from numpy.polynomial.hermite import poly2herm, herm2poly
    from numpy.polynomial.polynomial import polyder
    
    for _ in range(100): #1D
        h = np.random.rand(np.random.randint(1, 10))
        
        prediction = hermnpder(h)
        actual = poly2herm(polyder(herm2poly(h)))
        assert np.allclose(prediction, actual)
    
    for _ in range(100): #2D
        hs = np.random.rand(2, np.random.randint(1, 10))
        
        prediction = hermnpder(hs)
        actual = [poly2herm(polyder(herm2poly(h))) for h in hs]
        assert np.allclose(prediction, actual)

def test_hermnpderM():
    for _ in range(100):
        deg = np.random.randint(0, 10)
        h = np.random.rand(deg+1)
        
        prediction = hermnpderM(deg) @ h
        actual = hermnpder(h)
        assert np.allclose(prediction, actual)
