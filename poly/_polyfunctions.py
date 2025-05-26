from math import sumprod
from operator import neg, mul, pow, truediv
from functools import reduce
from itertools import repeat, count, chain, accumulate, islice
from vector import *
import sympy as sp
from sympy.abc import x as spx



__all__ = ['polyzero', 'polyone', 'polyx',
        'polyrand', 'polyrandn', 'polyfromroots',
        'polydeg', 'polyval', 'polycom',
        'polyaddc', 'polymul', 'polymulx', 'polypow', 'polydiv',
        'polyder', 'polyint',
        'polysympify', 'polyunsympify']



#creation stuff
#it's just veczero, but polyone and polyx are needed
#and therefore polyzero is also added again for consistency
polyzero = ()
"""Zero polynomial."""
polyone = (1, )
"""Constant one polynomial."""
polyx = (0, 1)
"""Identity ($f(x)=x$) polynomial."""

#For `polymono` use `vecbasis`

def polyrand(deg): #because deg+1
    return vecrand(deg+1)

def polyrandn(deg, normed=True, mu=0, sigma=1):
    return vecrandn(deg+1, normed=normed, mu=mu, sigma=sigma)

def polyfromroots(*rs):
    """Return the polynomial with the given roots."""
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    return polymul(*zip(map(neg, rs), repeat(1)))


#utility stuff
#For `polyeq` use `veceq`
#For `polytrim` use `vectrim`
#For `polyround` use `vecround`

def polydeg(p):
    """Return the degree of `p`.
    
    $\\deg(0)=-1$ is used for the empty zero polynomial
    (https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial).
    Doesn't handle leading zeros, use `vectrim` if needed.
    """
    return len(p) - 1


#evaluation
def polyval_naive(p, x):
    #Naive, N(N+1)/2 multiplications, N additions
    
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    #wouldn't need "if not p", but then the result may be of other type
    #return sumprod(p, map(pow, repeat(x), range(len(p))))
    
    #would work for iterables (without len)
    return sum(map(mul, p, map(pow, repeat(x), count())), start=type(x)(0))
    #return sum(map(mul, p, chain((type(x)(1),), map(pow, repeat(x), count(1)))), start=type(x)(0))

def polyval_iterative(p, x):
    #Iterative powers, 2N-1 multiplications, N additions
    return sum(map(mul, p, accumulate(repeat(x), mul, initial=type(x)(1))), start=type(x)(0))
    #return sum(map(mul, p, chain((type(x)(1),), accumulate(repeat(x), mul))), start=type(x)(0))

def polyval_horner(p, x):
    #Horner's method, N multiplications, N additions
    #https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency
    return reduce(lambda a, pi: a*x+pi, reversed(p), type(x)(0))

def polyval(p, x, method='horner'):
    """Return the value of `p` for the argument `x`."""
    match method:
        case 'naive':
            return polyval_naive(p, x)
        
        case 'iterative':
            return polyval_iterative(p, x)
        
        case 'horner':
            return polyval_horner(p, x)
        
        case _:
            raise ValueError

def polycom_naive(p, q):
    return vecadd(*(vecmul(pi, polypow(q, i)) for i, pi in enumerate(p)))

def polycom_iterative(p, q):
    return vecadd(*map(vecmul, p, accumulate(repeat(q), polymul, initial=polyone)))
    #return vecadd(*map(vecmul, p, chain((polyone,), accumulate(repeat(q), polymul))))

def polycom_horner(p, q):
    return reduce(lambda a, pi: polyaddc(polymul(a, q), pi), reversed(p), polyzero)

def polycom(p, q, method='horner'):
    """Return the composition $p(q)$."""
    match method:
        case 'naive':
            return polycom_naive(p, q)
        case 'iterative':
            return polycom_iterative(p, q)
        case 'horner':
            return polycom_horner(p, q)
        case _:
            return ValueError



#arithmetic stuff
#For `polyadd` use `vecadd`
#For `polysub` use `vecsub`
#For `polyscalarmul` use `vecmul`
#For `polyscalartruediv` use `vectruediv`
#For `polyscalarfloordiv` use `vecfloordiv`

def polyaddc(p, c):
    """Return the sum of a polynomial and a scalar."""
    it = iter(p)
    first = next(it, 0)
    return tuple(chain((first+c,), it))

def polymul_naive(p, q):
    r"""Return the product of two polynomials.
    
    Result: $\deg(pq)=\deg p+\deg q$, $\text{len}(pq)=\text{len}p+\text{len}q-1$.
    Operations
    - Additions: $(\text{len}p-1)(\text{len}q-1)=\deg p\deg q$
    - Multiplications: $\text{len}p\cdot\text{len}q=(\deg p+1)(\deg q+1)$
    The arguments must be sequences, not single exhaustible iterables.
    """
    p, q = tuple(p), tuple(q)
    if not p or not q:
        return polyzero
    r = [0] * (len(p)+len(q)-1)
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            r[i+j] += pi * qj
    return tuple(r)
    #r = [0] * (sum(map(len, ps)) - len(ps) + 1)
    #for i in product(*map(range, map(len, ps))):
    #    r[sum(i)] += prod(map(getitem, ps, i))
    #return tuple(r)

def _polymul_karatsuba(p, q):
    """Return $pq$ where $p$ & $q$ are of equal length and non empty."""
    if len(p) == 1:
        return (p[0] * q[0],)
    
    m = len(p) // 2
    pl, pu = p[:m], p[m:]
    ql, qu = q[:m], q[m:]
    
    z0 = _polymul_karatsuba(pl, ql)
    z2 = _polymul_karatsuba(pu, qu)
    z1 = _polymul_karatsuba(vecadd(pl, pu), vecadd(ql, qu))
    mid = vecsub(vecsub(z1, z0), z2)
    
    #return vecadd(z0, (0,)*m+mid, (0,)*(2*m)+z2)
    return z0[:m] + vecadd(z0[m:2*m], mid[:m]) + vecadd(z0[2*m:], mid[m:], z2)

def polymul_karatsuba(p, q):
    p, q = tuple(p), tuple(q)
    if not p or not q:
        return ()
    
    p, q = sorted((q, p), key=len) #p shorter, q longer
    ql, qu = q[:len(p)], q[len(p):]
    rl, ru = _polymul_karatsuba(p, ql), polymul_naive(p, qu)
    #return vecadd(rl, (0,)*len(p)+ru)
    return rl[:len(p)] + vecadd(rl[len(p):], ru)

def polymul(*ps, method='naive'):
    """Return the product of polynomials."""
    match method:
        case 'naive':
            return reduce(polymul_naive, ps, polyone)
        
        case 'karatsuba':
            return reduce(polymul_karatsuba, ps, polyone)
        
        case _:
            raise ValueError

def polymulx(p, n=1):
    """Return the product of a polynomial with its argument."""
    return tuple(chain((0,)*n, p))


def polypow_naive(p, n):
    return reduce(polymul, repeat(p, n), initial=polyone)

def polypow_binary(p, n):
    r = polyone
    while n:
        if n % 2 == 1:
            r = polymul(r, p)
        p = polymul(p, p)
        n //= 2
    return r

def polypow(p, n, method='binary'):
    """Return `p` raised to the `n`-th power."""
    match method:
        case 'naive':
            return polypow_naive(p, n)
        
        case 'binary':
            return polypow_binary(p, n)
        
        case _:
            raise ValueError


def polydiv(n, d):
    """Return $(q, r)$ such that $n=qd+r$."""
    #https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
    n, d, q = vectrim(n), vectrim(d), polyzero
    while n and polydeg(n)>=polydeg(d):
        t = vecbasis(len(n)-len(d), n[-1]/d[-1])
        q = vecadd(q, t)
        n = vectrim(vecsub(n, polymul(d, t))[:-1])
    return q, n

def polysqrt(p):
    p = tuple(p)
    if not p:
        return ()
    



#calculus stuff
def polyder(p, n=1):
    """Return the `n`-th derivative of `p`."""
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    for _ in range(n):
        p = vechadamard(count(1), islice(p, 1, None))
    return p

def polyint(p, n=1, c=0):
    """Return the `n`-th antiderivative of `p`.
    
    The integration constant `c` may be a sequence,
    or a scalar that is used in every integration step.
    """
    try:
        for i in range(n):
            p = (c[i],) + vechadamardtruediv(p, count(1))
        return p
    except TypeError:
        return polyint(p, n, (c,)*n)


#sympy
def polysympify(p, symbol=spx):
    return sp.Poly.from_list(tuple(reversed(p)), symbol)

def polyunsympify(p):
    return tuple(reversed(p.all_coeffs()))
