from math import sumprod
from operator import pow, mul, truediv, sub, eq
from functools import reduce
from itertools import zip_longest, repeat, starmap, chain



#constants
polyzero = ()
polyone = (1, )
polyx = (0, 1)

def mono(n, c=1):
    """Returns the monomial $cx^n$ as tuple."""
    return (0,)*(n-1) + (c,)


#utility
def polydeg(p):
    """Polynomial degree. Returns $\\deg(p)$,
    where $\\deg(0)=-1$ is used for the empty zero polynomial (https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial).
    Handles leading zeros"""
    return len(polytrim(p)) - 1

def polytrim(p, tol=0):
    """Polynomial trimming. Returns the polynomial, but with the leading
    coefficients, whos absolute values are less or equal to `tol`, removed."""
    while p and abs(p[-1])<=tol:
        p = p[:-1]
    return p

def polyeq(p, q):
    """Polynomial comparison. Returns `p==q`,
    ignoring leading zero coefficients."""
    return all(starmap(eq, zip_longest(p, q, fillvalue=0)))


#evaluation
def polyval(p, x):
    """Polynomial evaluation. Returns $p(x)$."""
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    return sumprod(p, map(pow, repeat(x), range(len(p)))) if p else type(x)(0)


#arithmetic
def polyadd(*ps):
    """Polynomial addition. Returns $\\sum_ip_i$. Handles empty sum."""
    return tuple(map(sum, zip_longest(*ps, fillvalue=0)))

def polysub(p, q):
    """Polynomial subtraction. Returns $p-q$."""
    return tuple(starmap(sub, zip_longest(p, q, fillvalue=0)))

def polymul(*ps):
    """Polynomial multiplication. Returns $\\prod_ip_i$.
    Handles empty product."""
    def _polymul(p, q):
        r = [0] * (len(p) + len(q) - 1)
        for i, pi in enumerate(p):
            for j, qj in enumerate(q):
                r[i+j] += pi * qj
        return r
    #wrap in tuple to be consistent with other functions
    #ternary is faster than initializer=polyone
    return tuple(reduce(_polymul, ps)) if ps else polyone
    #prolly faster for many arguments,
    #but isn't as clean and doesn't handle empty products:
    #r = [0] * (sum(len(pi) for pi in p) - len(p) + 1)
    #for i in product(*(range(len(pi)) for pi in p)):
    #    r[sum(i)] += prod(pi[ii] for pi, ii in zip(p, i))
    #return r

def polypow(p, n):
    """Polynomial exponentiation. Returns $p^n$. Handles `n=0`."""
    return polymul(*repeat(p, n)) #n=0 handled by empty product in polymul


#calculus
def polyder(p): #TODO: higher derivatives
    """Polynomial differentation. Returns $p'$."""
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    return tuple(map(mul, range(1, len(p)), p[1:]))

def polyint(p, c=0):
    """Polynomial integration. Returns $\\int p(x)dx$
    where `c` is the integration constant."""
    return tuple(chain((c,), map(truediv, p, range(1, len(p)+1))))



if __name__ == '__main__':
    import numpy as np
    from random import gauss, randint
    from tqdm import trange
    
    def gauss_tuple(n):
        return tuple(gauss() for _ in range(n))
    
    
    assert polydeg(()) == -1
    assert polydeg((0, )) == -1
    assert polydeg((0, 0)) == -1
    assert polydeg((1, )) == 0
    assert polydeg((1, 0)) == 0
    assert polydeg((1, 2)) == 1
    
    assert polytrim(polyzero) == polyzero
    assert polytrim((0, 0, 0)) == polyzero
    assert polytrim((1, 2, 3, 0, 0)) == (1, 2, 3)
    
    assert polyeq((), ()) == True
    assert polyeq((), (0,)) == True
    assert polyeq((), (1,)) == False
    assert polyeq((1,), (1,)) == True
    assert polyeq((1, 0), (1, )) == True
    assert polyeq((1, 0), (1, 2)) == False
    
    
    for _ in trange(10000, desc='polyval'):
        p, x = gauss_tuple(randint(1, 10)), gauss()
        assert np.isclose(polyval(p, x),
                np.polynomial.polynomial.polyval(x, p))
    assert polyval(polyzero, 1) == 0
    assert type(polyval(polyzero, 1.0)) == float
    
    
    for _ in trange(10000, desc='polyadd'):
        ps = [gauss_tuple(randint(1, 10)) for _ in range(randint(1, 4))]
        assert np.allclose(polyadd(*ps),
                reduce(np.polynomial.polynomial.polyadd, ps))
    assert polyeq(polyadd(), polyzero)
    assert polyeq(polyadd(p), p)
    assert polyeq(polyadd(p, polyzero), p)
    
    for _ in trange(10000, desc='polysub'):
        p, q = gauss_tuple(randint(1, 10)), gauss_tuple(randint(1, 10))
        assert np.allclose(polysub(p, q),
                np.polynomial.polynomial.polysub(p, q))
    assert polyeq(polysub(p, polyzero), p)
    assert polyeq(polysub(polyzero, p), tuple(-pi for pi in p))
    
    for _ in trange(10000, desc='polymul'):
        ps = [gauss_tuple(randint(1, 10)) for _ in range(randint(1, 4))]
        assert np.allclose(polymul(*ps),
                reduce(np.polynomial.polynomial.polymul, ps))
    assert polyeq(polymul(), polyone)
    assert polyeq(polymul(p), p)
    assert polyeq(polymul(p, polyzero), polyzero)
    assert polyeq(polymul(p, polyone), p)
    
    for _ in trange(1000, desc='polypow'):
        p, n = gauss_tuple(randint(1, 10)), randint(0, 20)
        assert np.allclose(polypow(p, n),
                np.polynomial.polynomial.polypow(p, n))
    assert polyeq(polypow(p, 0), polyone)
    assert polyeq(polypow(p, 1), p)
    
    
    for _ in trange(1000, desc='polyder'):
        p = gauss_tuple(randint(1, 10))
        assert np.allclose(polyder(p),
                np.polynomial.polynomial.polyder(p))
    assert polyeq(polyder(polyzero), polyzero)
    assert polyeq(polyder(polyone), polyzero)
    
    for _ in trange(1000, desc='polyint'):
        p, c = gauss_tuple(randint(1, 10)), gauss()
        assert np.allclose(polyint(p, c),
                np.polynomial.polynomial.polyint(p, k=c))
    assert polyeq(polyint(polyzero), polyzero)
    assert polyeq(polyint(polyzero, 1), polyone)
    assert polyeq(polyint(polyone), polyx)
