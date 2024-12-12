import numpy as np
from functools import reduce
from itertools import repeat
from math import sumprod
from vector import vecnpzero, vecnpbasis, vecnptrim, vecnpeq, vecnpadd



__all__ = ['polynpzero', 'polynpone', 'polynpmono', 'polynpfromroots',
        'polynpval', 'polynpcom', 'polynpmul', 'polynppow',
        'polynpder', 'polynpmder', 'polynpint', 'polynpmint']



#creation stuff
def polynpzero(d=None):
    """Return `d` zero polynomials.
    
    The retured value is a `(d, 1)`-array of zeros if `d` is not `None`
    or `[0]` otherwise.
    """
    return vecnpzero(d)

def polynpone(d=None):
    """Return `d` constant polynomials.
    
    The retured value is a `(d, 1)`-array of ones if `d` is not `None`
    or `[0]` otherwise.
    """
    #same dtype as numpy.polynomial.polynomial.polyone
    return np.ones(1 if d is None else (d, 1), dtype=np.int64)

def polynpmono(n, c=1, d=None):
    """Return `d` many $cx^n$ polynomials (`n*[0]+[c]`).
    
    The retured value is a `(d, i+1)`-array if `d` is not `None`
    or `(i+1,)` otherwise.
    """
    return vecnpbasis(n, c=c, d=d)


def polynpfromroots(r):
    r = np.asarray(r)
    return polynpmul(*(
            np.pad(-ri.T, ((0, 0),)*(r.ndim-1)+((0, 1),), constant_values=1)
            for ri in r.T[:,np.newaxis]))


#utility
#polydeg not implemented as numpy uses [0] as polyzero and not []
#it therefore deviates from the clean approach in poly.py
#and polydeg is also not implemented in numpy.polynomial.polynomial


#evaluation
def polynpval(p, x):
    """Return the value of `p` for the argument `x`."""
    p = np.asarray(p)
    return sumprod(p, map(pow, repeat(x), range(p.shape[-1])))

def polynpcom(p, q):
    """Polynomial composition. Returns $p(q)$."""
    t, r = np.ones((p.shape[0], 1)), np.zeros((p.shape[0], 1))
    for pi in p.T:
        r = vecnpadd(r, t*pi[:,np.newaxis])
        t = polynpmul(t, q)
    return r


#arithmetic
def _polynpmul(p, q):
    r = [0] * (p.shape[-1] + q.shape[-1] - 1)
    for i, pi in enumerate(p.T):
        for j, qj in enumerate(q.T):
            r[i+j] += pi * qj
    return np.array(r).T

def polynpmul(*ps):
    if not ps:
        return polynpzero()
    return reduce(_polynpmul, map(np.asarray, ps))

def polynppow(p, n):
    p = np.asarray(p)
    if not n:
        return polynpone(p.shape[0] if p.ndim==2 else None)
    return polynpmul(*repeat(p, n))


#calculus
def polynpder(p, n=1):
    p = np.asarray(p)
    for _ in range(n):
        p = p[...,1:] * np.arange(1, p.shape[-1])
    return p

def polynpmder(deg):
    D = np.zeros((deg, deg+1), dtype=object)
    D.flat[1::D.shape[1]+1] = range(1, deg+1)
    return D

def polynpint(p, n=1, c=0):
    p = np.asarray(p)
    try:
        for i in range(n):
            p = np.pad(p/np.arange(1, p.shape[-1]+1),
                    (1, 0) if p.ndim==1 else ((0, 0), (1, 0)),
                    constant_values=c[i])
        return p
    except TypeError:
        return polynpint(p, n, (c,)*n)

def polynpmint(deg):
    I = np.zeros((deg+2, deg+1), dtype=np.float64)
    I.flat[I.shape[1]::I.shape[1]+1] = 1 / np.arange(1, deg+2)
    return I
