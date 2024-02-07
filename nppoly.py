import numpy as np
from functools import reduce
from itertools import repeat



#constants
def polyzero(d):
    """Returns `d` many zero polynomials (`[0]`)."""
    return np.zeros((d, 1))

def polyone(d):
    """Returns `d` many one polynomials (`[1]`)."""
    return np.ones((d, 1))

def polymono(d, n):
    """Returns `d` many $x^n$ polynomials (`n*[0]+[1]`)."""
    return np.pad(polyone(d), ((0, 0), (n, 0)))


#utility
#polydeg not implemented as numpy uses [0] as polyzero and not []
#it therefore deviates from the clean approach in poly.py
#and polydeg is also not implemented in numpy.polynomial.polynomial

def polytrim(p, tol=0):
    """Polynomial trimming. Returns the polynomial, but with the leading
    coefficients, whos absolute values are less or equal to `tol`, removed."""
    while p.shape[1]>1 and any(abs(p[:,-1])<=tol):
        p = p[:,:-1]
    return p

def polyeq(p, q):
    """Polynomial comparison. Returns `p==q`,
    ignoring leading zero coefficients."""
    return np.array_equal(polytrim(p), polytrim(q))


#evaluation
#numpy.polynomial.polynomial.polyval already vectorized

def polycom(p, q):
    """Polynomial composition. Returns $p(q)$."""
    t, r = np.ones((p.shape[0], 1)), np.zeros((p.shape[0], 1))
    for pi in p.T:
        r = polyadd(r, t*pi[:,np.newaxis])
        t = polymul(t, q)
    return r


#arithmetic
def polyadd(*ps):
    """Polynomial addition. Returns $\\sum_ip_i$. Doesn't handle empty sum."""
    #r = ps[0] padded and then adding ps[1:] is slower
    r = np.zeros((ps[0].shape[0], max(p.shape[1] for p in ps)))
    for p in ps:
        r[:,:p.shape[1]] += p
    return r

def polysub(p, q):
    """Polynomial subtraction. Returns $p-q$."""
    #r = p padded and then r[:,:q.shape[1]]-=q is slower
    r = np.zeros((p.shape[0], max(p.shape[1], q.shape[1])))
    r[:,:p.shape[1]] += p
    r[:,:q.shape[1]] -= q
    return r

def polymul(*p):
    """Polynomial multiplication. Returns $\\prod_ip_i$.
    Doesn't handle empty product."""
    def _polymul(p, q):
        r = np.zeros((p.shape[0], p.shape[1]+q.shape[1]-1))
        for i, pi in enumerate(p.T):
            for j, qj in enumerate(q.T):
                r[:,i+j] += pi * qj
        return r
    return reduce(_polymul, p)

def polypow(p, n):
    """Polynomial exponentiation. Returns $p^n$. Handles `n=0`."""
    return polymul(*repeat(p, n)) if n!=0 else polyone(p.shape[0])


#calculus
def polyder(p):
    """Polynomial differentation. Returns $p'$."""
    return p[:,1:] * np.arange(1, p.shape[1])

def polyint(p, c=0):
    """Polynomial integration. Returns $\\int p(x)dx$
    where `c` is the integration constant."""
    return np.pad(p/np.arange(1, p.shape[1]+1), ((0, 0), (1, 0)),
            constant_values=c)



if __name__ == '__main__':
    from tqdm import trange
    
    
    assert np.array_equal(polyzero(1), np.array([[0]]))
    assert np.array_equal(polyzero(2), np.array([[0], [0]]))
    assert np.array_equal(polyzero(3), np.array([[0], [0], [0]]))
    
    assert np.array_equal(polyone(1), np.array([[1]]))
    assert np.array_equal(polyone(2), np.array([[1], [1]]))
    assert np.array_equal(polyone(3), np.array([[1], [1], [1]]))
    
    assert np.array_equal(polyone(1), np.array([[1]]))
    assert np.array_equal(polyone(2), np.array([[1], [1]]))
    assert np.array_equal(polyone(3), np.array([[1], [1], [1]]))
    
    assert np.array_equal(polymono(1, 0), np.array([[1]]))
    assert np.array_equal(polymono(2, 0), np.array([[1], [1]]))
    assert np.array_equal(polymono(1, 1), np.array([[0, 1]]))
    assert np.array_equal(polymono(2, 1), np.array([[0, 1], [0, 1]]))
    
    
    assert np.array_equal(polytrim(np.array([[0]])), np.array([[0]]))
    assert np.array_equal(polytrim(np.array([[1, 0], [2, 0]])), np.array([[1], [2]]))
    
    assert polyeq(np.array([[0]]), np.array([[0]])) == True
    assert polyeq(np.array([[0]]), np.array([[0, 0]])) == True
    assert polyeq(np.array([[0]]), np.array([[1]])) == False
    assert polyeq(np.array([[1]]), np.array([[1]])) == True
    assert polyeq(np.array([[1, 0]]), np.array([[1]])) == True
    assert polyeq(np.array([[1, 0]]), np.array([[1, 2]])) == False
    
    
    for _ in trange(1000, desc='polycom'):
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        q = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polycom(p, q)
        actual = [np.poly1d(pi[::-1])(np.poly1d(qi[::-1])).c[::-1] for pi, qi in zip(p, q)]
        assert np.allclose(prediction, actual)
    
    
    for _ in trange(1000, desc='polyadd'):
        n = np.random.randint(1, 4)
        ps = [np.random.normal(size=(20, np.random.randint(1, 10))) for _ in range(n)]
        prediction = polyadd(*ps)
        actual = [reduce(np.polynomial.polynomial.polyadd, rows) for rows in zip(*ps)]
        assert np.allclose(prediction, actual)
    
    for _ in trange(1000, desc='polysub'):
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        q = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polysub(p, q)
        actual = [np.polynomial.polynomial.polysub(pi, qi) for pi, qi in zip(p, q)]
        assert np.allclose(prediction, actual)
    
    for _ in trange(1000, desc='polymul'):
        n = np.random.randint(1, 4)
        ps = [np.random.normal(size=(20, np.random.randint(1, 10))) for _ in range(n)]
        prediction = polymul(*ps)
        actual = [reduce(np.polynomial.polynomial.polymul, rows) for rows in zip(*ps)]
        assert np.allclose(prediction, actual)
    
    for _ in trange(1000, desc='polypow'):
        n = np.random.randint(0, 5)
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polypow(p, n)
        actual = np.apply_along_axis(
                lambda pi: np.polynomial.polynomial.polypow(pi, n), 1, p)
        assert np.allclose(prediction, actual)
    
    
    for _ in trange(1000, desc='polyder'):
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polyder(p)
        actual = np.apply_along_axis(
                np.polynomial.polynomial.polyder, 1, p)
        assert np.allclose(prediction, actual)
    
    for _ in trange(1000, desc='polyint'):
        c = np.random.normal()
        p = np.random.normal(size=(20, np.random.randint(1, 10)))
        prediction = polyint(p, c)
        actual = np.apply_along_axis(
                lambda p: np.polynomial.polynomial.polyint(p, k=c), 1, p)
        assert np.allclose(prediction, actual)
