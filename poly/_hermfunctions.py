from math import sumprod, factorial, comb
from operator import neg, mul, pow, truediv
from fractions import Fraction
from functools import reduce, cache
from itertools import repeat, chain, accumulate, islice, count, tee, zip_longest
from vector import vecbasis, vecrand, vecadd, vecsub, vecmul, vectrim, veceq, vechadamard
import sympy as sp
from sympy.abc import x as spx

from ._polyfunctions import polyzero, polyone, polyval, polyaddc, polymulx



__all__ = ['herm', 'herm2poly', 'hermmono', 'poly2herm',
           'hermval',
           'hermdot', 'hermspdot',
           'hermmulx', 'hermmul', 'hermpow', 'hermmulpow',
           'hermder', 'hermint',
           'hermsympify']



#@cache
def herm_recursive(n):
    #https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation
    if n == 0:
        return (1,)
    elif n == 1:
        return (0, 2)
    else:
        return vecmul(2, vecsub(polymulx(herm_recursive(n-1)),
                                vecmul(n-1, herm_recursive(n-2))))

def herm_iterative(n):
    if n == 0:
        return (1,)
    Him1, Hi = (1,), (0, 2)
    for i in range(2, n+1):
        Him1, Hi = Hi, vecmul(2, vecsub(polymulx(Hi), vecmul(i-1, Him1)))
    return Hi

def herm_explicit(n):
    r = [0] * (n+1)
    r[-1] = 2**n
    for m in range(1, n//2+1):
        r[n-2*m] = -r[n-2*(m-1)] * (n-2*m+1)*(n-2*m+2) // (4*m)
    return tuple(r)

def herm(n, method='iterative'):
    """Return the `n`-th Hermite polynomial."""
    match method:
        case 'recursive':
            return herm_recursive(n)
        
        case 'iterative':
            return herm_iterative(n)
        
        case 'explicit':
            return herm_explicit(n)
        
        case _:
            raise ValueError


#conversion stuff
def herm2poly_naive(h):
    return vecadd(*(vecmul(hi, herm(i)) for i, hi in enumerate(h)))

def herm2poly_clenshaw(h):
    h = tuple(h)
    if not h:
        return polyzero
    a, b = polyzero, polyzero
    for n, hn in reversed(tuple(enumerate(h[1:], 1))):
        a, b = polyaddc(vecsub(vecmul(2, polymulx(a)), vecmul(2*(n+1), b)), hn), a
    return polyaddc(vecsub(vecmul(2, polymulx(a)), vecmul(2, b)), h[0])

def herm2poly(h, method='naive'):
    match method:
        case 'naive':
            return herm2poly_naive(h)
        
        case 'clenshaw':
            return herm2poly_clenshaw(h)
        
        case _:
            raise ValueError

def hermmono_recursive(n):
    if n == 0:
        return (Fraction(1),)
    return hermmulx(hermmono(n-1))

def hermmono_iterative(n):
    h = (Fraction(1),)
    for _ in range(n):
        h = hermmulx(h)
    return h

def hermmono_explicit(n):
    r = [Fraction()] * (n + 1)
    #for m in range(n//2+1):
    #    r[n-2*m] = Fraction(factorial(n), 2**n * factorial(m)*factorial(n-2*m))
    r[n] = Fraction(1, 2**n)
    for m in range(1, n//2+1):
        r[n-2*m] = r[n-2*(m-1)] * ((n-2*m+1)*(n-2*m+2)) / m
    return r

def hermmono(n, method='explicit'):
    match method:
        case 'recursive':
            return hermmono_recursive(n)
        
        case 'iterative':
            return hermmono_iterative(n)
        
        case 'explicit':
            return hermmono_explicit(n)
        
        case _:
            raise ValueError

def poly2herm_naive(p):
    return vecadd(*(vecmul(pi, hermmono(i)) for i, pi in enumerate(p)))

def poly2herm(p):
    #TODO: h = sum_i int_\mathbb{R}exp(-x^2)p(x)H_i(x)dx H_i
    #int_\mathbb{R}exp(-x^2)p(x)H_i(x)dx = sum_n p_n int_\mathbb{R}exp(-x^2)x^nH_i(x)dx
    h = [0]*len(p)
    for i in reversed(range(len(p))):
        Hi = herm(i)
        t = vecmul(p[i]/Hi[-1], Hi)
        h[i] += t[i]
        p = vecsub(p, t)
    return tuple(h)


#evaluation stuff
def hermval_naive(h, x):
    return polyval(herm2poly(h), x)

def hermval_clenshaw(h, x):
    h = tuple(h)
    if not h:
        return type(x)(0)
    a, b = 0, 0
    for n, hn in reversed(tuple(enumerate(h[1:], 1))):
        a, b = hn + 2*x*a - 2*(n+1)*b, a
    return h[0] + 2*x*a - 2*b

def hermval(h, x, method='clenshaw'):
    """Return the value of `h` for the argument `x`."""
    match method:
        case 'naive':
            return hermval_naive(h, x)
        
        case 'clenshaw':
            return hermval_clenshaw(h, x)
        
        case _:
            raise ValueError


def hermdot(g, h):
    """Return $\frac{1}{\sqrt{\pi}}\int_\mathbb{R}e^{-x^2}g(x)h(x)dx$."""
    return sum(gi*hi * 2**i*factorial(i) for i, (gi, hi) in enumerate(zip(g, h)))

def hermspdot(g, h):
    """Return $\int_\mathbb{R}e^{-x^2}g(x)h(x)dx$."""
    return sp.sqrt(sp.pi) * hermdot(g, h)


def hermmulx(h):
    """Return $xh(x)$."""
    a, b = tee(h)
    _, c = next(a, None), next(a, 0)
    return tuple(chain((c,), (hkm1/2+(k+1)*hkp1 for k, (hkp1, hkm1) in
            enumerate(zip_longest(a, b, fillvalue=0), 1))))

def _hermmul(a, b):
    """Return $ab$."""
    a, b = tuple(a), tuple(b)
    if not a or not b:
        return polyzero
    r = [0] * (len(a)+len(b)-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            for k in range(min(i, j)+1):
                r[i+j-2*k] += 2**k*factorial(k)*comb(i, k)*comb(j, k) * ai*bj
    
    return tuple(r)

def hermmul(*hs):
    return reduce(_hermmul, hs, polyone)

def hermpow_naive(h, n):
    """Return `h` raised to the `n`-th power."""
    return reduce(polymul, repeat(p, n), polyone)

def hermpow_binary(h, n):
    """Return `h` raised to the `n`-th power."""
    r = polyone
    while n:
        if n % 2 == 1:
            r = hermmul(r, h)
        h = hermmul(h, h)
        n //= 2
    return r

def hermpow(h, n, method='binary'):
    match method:
        case 'naive':
            return hermpow_naive(h, n)
        case 'binary':
            return hermpow_binary(h, n)
        case _:
            raise ValueError

def hermmulpow(alpha):
    """Return $H^\alpha=\prod_iH_i^{\alpha_i}$."""
    r = polyone
    for i, ai in enumerate(alpha):
        r = hermmul(hermpow(vecbasis(i), ai))
    return r


#calculus stuff
def hermder(h, n=1):
    """Return the `n`-th derivative of `h`."""
    for _ in range(n):
        h = vechadamard(islice(h, 1, None), count(2, 2))
    return h

def hermint(h, n=1):
    """Return the `n`-th antiderivative of `h`."""
    for _ in range(n):
        h = (0,) + vechadamardtruediv(h, count(2, 2))
    return h


def hermsympify(h, symbol=spx):
    return sum((hi*sp.hermite(i, symbol) for i, hi in enumerate(h)), start=sp.Integer(0))
