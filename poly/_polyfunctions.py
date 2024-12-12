from math import sumprod
from operator import neg, mul, pow, truediv
from functools import reduce
from itertools import repeat, chain, accumulate
from vector import vecbasis, vecrand, vecadd, vecsub, vecmul, vectrim, veceq



__all__ = ['polyzero', 'polyone', 'polyx',
        'polyrand', 'polyrandn', 'polyfromroots',
        'polydeg', 'polyval', 'polycom',
        'polymul', 'polymulx', 'polypow', 'polydiv',
        'polyder', 'polyint']



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

def polyrandn(deg, normed=True):
    return vecrandn(deg+1, normed=normed)

def polyfromroots(rs):
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
def polyval(p, x, method='horner'):
    """Return the value of `p` for the argument `x`."""
    if not p:
        return type(x)(0)
    
    match method:
        case 'naive':
            #Naive, N(N+1)/2 multiplications, N additions
            #https://docs.python.org/3/library/itertools.html#itertools-recipes
            return sumprod(p, map(pow, repeat(x), range(len(p))))
        
        case 'iterative':
            #Iterative powers, 2N-1 multiplications, N additions
            return sumprod(p, accumulate(repeat(x, len(p)-1), mul, initial=1))
        
        case 'horner':
            #Horner's method, N multiplications, N additions
            #https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency
            return reduce(lambda a, c: a*x+c, reversed(p))

def polycom(p, q):
    """Return the composition $p(q)$."""
    if not p:
        return polyzero
    #iterative, because naive would use polypow, that itself is iterative
    qei, r = polyone, polyzero
    for pi in p[:-1]:
        r = vecadd(r, vecmul(pi, qei))
        qei = polymul(qei, q) #would be unnecessary in last iteration
    return vecadd(r, vecmul(p[-1], qei))


#arithmetic stuff
#For `polyadd` use `vecadd`
#For `polysub` use `vecsub`
#For `polyscalarmul` use `vecmul`
#For `polyscalartruediv` use `vectruediv`
#For `polyscalarfloordiv` use `vecfloordiv`

def _polymul(p, q):
    r = [0]*(len(p)+len(q)-1) if p and q else ()
    for i, pi in enumerate(p):
        for j, qj in enumerate(q):
            r[i+j] += pi * qj
    return r

def polymul(*ps):
    """Return the product of polynomials."""
    return tuple(reduce(_polymul, ps, polyone))
    #r = [0] * (sum(map(len, ps)) - len(ps) + 1)
    #for i in product(*map(range, map(len, ps))):
    #    r[sum(i)] += prod(map(getitem, ps, i))
    #return tuple(r)
    #TODO: Karatsuba

def polymulx(p):
    return tuple(chain((0,), p))

def polypow(p, n):
    """Return `p` raised to the `n`-th power."""
    return polymul(*repeat(p, n)) #n=0 handled by empty product in polymul

def polydiv(n, d):
    """Return $(q, r)$ such that $n=qd+r$."""
    #https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
    n, d, q = vectrim(n), vectrim(d), polyzero
    while n and polydeg(n)>=polydeg(d):
        t = vecbasis(len(n)-len(d), n[-1]/d[-1])
        q = vecadd(q, t)
        n = vectrim(vecsub(n, polymul(d, t))[:-1])
    return q, n


#calculus stuff
def polyder(p, n=1):
    """Return the `n`-th derivative of `p`."""
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    for _ in range(n):
        p = tuple(map(mul, range(1, len(p)), p[1:]))
    return p

def polyint(p, n=1, c=0):
    """Return the `n`-th antiderivative of `p`.
    
    The integration constant `c` may be a sequence,
    or a scalar that is used in every integration step.
    """
    try:
        for i in range(n):
            p = tuple(chain((c[i],), map(truediv, p, range(1, len(p)+1))))
        return p
    except TypeError:
        return polyint(p, n, (c,)*n)
