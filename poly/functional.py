from operator import neg, mul, pow
from functools import reduce
from itertools import repeat, count, islice
from vector import *
import sympy as sp
from sympy.abc import x as spx



__all__ = (#creation
           'polyzero', 'polyone', 'polyx', 'polymono',
           'polyrand', 'polyrandn', 'polyfromroots',
           #utility
           'polyeq', 'polytrim', 'polyround', 'polydeg',
           #evaluation
           'polyval', 'polyval_naive', 'polyval_iterative', 'polyval_horner', 'polyvalgen', 'polyvalzero',
           'polycom', 'polycom_naive', 'polycom_iterative', 'polycom_horner', 'polycomgen',
           'polyshift',
           #arithmetic
           'polypos', 'polyneg',
           'polyadd', 'polyaddc', 'polysub', 'polyscalarmul',
           'polyscalartruediv', 'polyscalarfloordiv',
           'polyscalarmod', 'polyscalardivmod',
           'polymul', 'polymul_naive', 'polymul_karatsuba',
           'polymulx',
           'polypow', 'polypow_naive', 'polypow_binary',
           #calculus
           'polyder', 'polyantider',
           #sympy
           'polysympify', 'polyunsympify')



#creation
polyzero = ()
"""Zero polynomial.

$$
    0
$$

An empty tuple: `()`.

Notes
-----
Why give the zero polynomial a distinguished representation (not just `[0]` like [`numpy.polynomial`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyzero.html))?

The additive neutral element seems worth handling exceptionally.
It is mathematically different (different degree)
and result in functions like [`polymul`][poly.functional.polymul] being more time and memory efficient.

See also
--------
- for any degree: [`polymono`][poly.functional.polymono]
- wraps: [`vector.veczero`](https://goessl.github.io/vector/functional/#vector.functional.veczero)
"""
polyone = (1, )
"""Constant one polynomial.

$$
    1
$$

A tuple with a single one: `(1,)`.

See also
--------
- for any degree: [`polymono`][poly.functional.polymono]
"""
polyx = (0, 1)
"""Identity polynomial.

$$
    x
$$

A tuple with a zero and a one `(0, 1)`.

See also
--------
- for any degree: [`polymono`][poly.functional.polymono]
"""

def polymono(n, c=1):
    """Return a monomial of degree `n`.
    
    $$
        cx^n
    $$
    
    Returns a tuple with `n` zeros followed by `c` or [`polyzero`][poly.functional.polyzero] if $n<0$.
    
    See also
    --------
    - constant for $n=-1$: [`polyzero`][poly.functional.polyzero]
    - constant for $n=0$: [`polyone`][poly.functional.polyone]
    - constant for $n=1$: [`polyx`][poly.functional.polyx]
    - wraps: [`vector.vecbasis`](https://goessl.github.io/vector/functional/#vector.functional.vecbasis)
    """
    return vecbasis(n, c) if n>=0 else polyzero

def polyrand(n):
    r"""Return a random polynomial of degree `n`.
    
    $$
        \sum_{k=0}^na_kx^k \quad a_k\sim\mathcal{U}([0, 1[)
    $$
    
    The coefficients are sampled from a uniform distribution in `[0, 1[`.
    
    See also
    --------
    - wraps: [`vector.vecrand`](https://goessl.github.io/vector/functional/#vector.functional.vecrand)
    """
    return vecrand(n+1)

def polyrandn(n, normed=True, mu=0, sigma=1):
    r"""Return a random polynomial of degree `n`.
    
    $$
        \sum_{k=0}^na_kx^k \quad a_k\sim\mathcal{N}(\mu, \sigma)
    $$
    
    The coefficients are sampled from a normal distribution.
    
    See also
    --------
    - wraps: [`vector.vecrandn`](https://goessl.github.io/vector/functional/#vector.functional.vecrandn)
    """
    return vecrandn(n+1, normed=normed, mu=mu, sigma=sigma)

def polyfromroots(*xs):
    r"""Return the polynomial with the given roots.
    
    $$
        \prod_k(x-x_k)
    $$
    
    References
    ----------
    - <https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.polynomial_from_roots>
    """
    return polymul(*zip(map(neg, xs), repeat(1)))


#utility
def polyeq(p, q):
    r"""Return if two polynomials are equal.
    
    $$
        p \overset{?}{=} q
    $$
    
    See also
    --------
    - wraps: [`vector.veceq`](https://goessl.github.io/vector/functional/#vector.functional.veceq)
    """
    return veceq(p, q)

def polytrim(p, tol=1e-9):
    """Remove all leading near zero (`abs(p_i)<=tol`) coefficients.
    
    See also
    --------
    - wraps: [`vector.vectrim`](https://goessl.github.io/vector/functional/#vector.functional.vectrim)
    """
    return vectrim(p, tol)

def polyround(p, ndigits=None):
    r"""Round all coefficients to the given precision.
    
    $$
        \sum_k\text{round}_\text{ndigits}(a_k)x^k
    $$
    
    See also
    --------
    - wraps: [`vector.vecround`](https://goessl.github.io/vector/functional/#vector.functional.vecround)
    """
    return vecround(p, ndigits)

def polydeg(p):
    r"""Return the degree of a polynomial.
    
    $$
        \deg(p)
    $$
    
    $\deg(0)=-1$ is used for the empty zero polynomial.
    
    Doesn't handle leading zeros, use [`polytrim`][poly.functional.polytrim]
    if needed.
    
    `p` must support `len(p)`.
    
    Notes
    -----
    $\deg(0)=-\infty$ is more commonly used but the expected return type is an
    `int` and `-math.inf` is of type `float`. Therefore the $\deg(0)=-1$
    convention was choosen to conserve the return type.
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial>
    """
    return len(p) - 1


#evaluation
def polyval_naive(p, x):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Uses naive repeated multiplication to calculate monomials individually.
    
    See also
    --------
    - for any implementation: [`polyval`][poly.functional.polyval]
    - other implementations:
    [`polyval_iterative`][poly.functional.polyval_iterative],
    [`polyval_horner`][poly.functional.polyval_horner]
    - for polynomial arguments: [`polycom_naive`][poly.functional.polycom_naive]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency>
    - <https://en.wikipedia.org/wiki/Polynomial_evaluation>
    """
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    #wouldn't need "if not p", but then the result may be of other type
    #return sumprod(p, map(pow, repeat(x), range(len(p))))
    
    #works for iterables (without len) and is type-safe
    return sum(map(mul, p, map(pow, repeat(x), count())), start=type(x)(0))

def polyvalgen(x):
    r"""Yield the powers of the value `x`.
    
    $$
        (x^n)_{n\in\mathbb{N}_0} = (1, x, x^2, \dots)
    $$
    
    Uses iterative multiplication to calculate monomials consecutively.
    
    See also
    --------
    - used in: [`polyval_iterative`][poly.functional.polyval_iterative]
    - for polynomial arguments: [`polycomgen`][poly.functional.polycomgen]
    """
    yield (y := type(x)(1))
    while True:
        y *= x
        yield y

def polyval_iterative(p, x):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Uses iterative multiplication to calculate monomials consecutively.
    
    See also
    --------
    - for any implementation: [`polyval`][poly.functional.polyval]
    - other implementations:
    [`polyval_naive`][poly.functional.polyval_naive],
    [`polyval_horner`][poly.functional.polyval_horner]
    - uses: [`polyvalgen`][poly.functional.polyvalgen]
    - for polynomial arguments: [`polycom_iterative`][poly.functional.polycom_iterative]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency>
    """
    return sum(map(mul, p, polyvalgen(x)), start=type(x)(0))

def polyval_horner(p, x):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Uses Horner's method.
    
    `p` must be reversible
    
    See also
    --------
    - for any implementation: [`polyval`][poly.functional.polyval]
    - other implementations:
    [`polyval_naive`][poly.functional.polyval_naive],
    [`polyval_iterative`][poly.functional.polyval_iterative]
    - for polynomial arguments: [`polycom_horner`][poly.functional.polycom_horner]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Horner%27s_method>
    """
    return reduce(lambda a, pi: a*x+pi, reversed(p), type(x)(0))

def polyval(p, x, method='horner'):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Available methods are
    
    - [`naive`][poly.functional.polyval_naive],
    - [`iterative`][poly.functional.polyval_iterative] &
    - [`horner`][poly.functional.polyval_horner] (`p` must be reversible).
    
    See also
    --------
    - implementations: [`polyval_naive`][poly.functional.polyval_naive],
    [`polyval_iterative`][poly.functional.polyval_iterative],
    [`polyval_horner`][poly.functional.polyval_horner]
    - for $x=0$: [`polyvalzero`][poly.functional.polyvalzero]
    - for consecutive monomials: [`polyvalgen`][poly.functional.polyvalgen]
    - for polynomial arguments: [`polycom`][poly.functional.polycom]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency>
    - <https://en.wikipedia.org/wiki/Polynomial_evaluation>
    """
    match method:
        case 'naive':
            return polyval_naive(p, x)
        case 'iterative':
            return polyval_iterative(p, x)
        case 'horner':
            return polyval_horner(p, x)
        case _:
            raise ValueError('Invalid method')

def polyvalzero(p):
    """Return the value of polynomial `p` evaluated at point 0.
    
    $$
        p(0)
    $$
    
    More efficient than `polyval(p, 0)`.
    
    See also
    --------
    - for any argument: [`polyval`][poly.functional.polyval]
    """
    return next(iter(p), 0)

def polycom_naive(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q \qquad \text{or similarly} \qquad p(q)
    $$
    
    Uses naive repeated multiplication to calculate monomials individually.
    
    `q` must be a sequence.
    
    See also
    --------
    - for any implementation: [`polycom`][poly.functional.polycom]
    - other implementations:
    [`polycom_iterative`][poly.functional.polycom_iterative],
    [`polycom_horner`][poly.functional.polycom_horner]
    - based on: [`polyval_naive`][poly.functional.polyval_naive]
    """
    return polyadd(*(polyscalarmul(pi, polypow(q, i)) for i, pi in enumerate(p)))

def polycomgen(p):
    r"""Yield the powers of the polynomial `p`.
    
    $$
        (p^n)_{n\in\mathbb{N}_0} = (1, p, p^2, \dots)
    $$
    
    Uses iterative multiplication to calculate powers consecutively.
    
    `p` must be a sequence.
    
    See also
    --------
    - used in: [`polycom_iterative`][poly.functional.polycom_iterative]
    - based on: [`polyvalgen`][poly.functional.polyvalgen]
    """
    yield (q := polyone)
    while True:
        q = polymul(q, p)
        yield q

def polycom_iterative(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q \qquad \text{or similarly} \qquad p(q)
    $$
    
    Uses iterative multiplication to calculate monomials consecutively.
    
    `q` must be a sequence.
    
    See also
    --------
    - for any implementation: [`polycom`][poly.functional.polycom]
    - other implementations:
    [`polycom_naive`][poly.functional.polycom_naive],
    [`polycom_horner`][poly.functional.polycom_horner]
    - uses: [`polycomgen`][poly.functional.polycomgen]
    - based on: [`polyval_iterative`][poly.functional.polyval_iterative]
    """
    return polyadd(*map(polyscalarmul, p, polycomgen(q)))

def polycom_horner(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q \qquad \text{or similarly} \qquad p(q)
    $$
    
    Uses Horner's method.
    
    `q` must be reversible.
    
    See also
    --------
    - for any implementation: [`polycom`][poly.functional.polycom]
    - other implementations:
    [`polycom_naive`][poly.functional.polycom_naive],
    [`polycom_iterative`][poly.functional.polycom_iterative]
    - based on: [`polyval_horner`][poly.functional.polyval_horner]
    """
    return reduce(lambda a, pi: polyaddc(polymul(a, q), pi), reversed(p), polyzero)

def polycom(p, q, method='iterative'):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q \qquad \text{or similarly} \qquad p(q)
    $$
    
    Available methods are
    
    - [`naive`][poly.functional.polycom_naive],
    - [`iterative`][poly.functional.polycom_iterative] &
    - [`horner`][poly.functional.polycom_horner] (`p` must be reversible).
    
    See also
    --------
    - implementations: [`polycom_naive`][poly.functional.polycom_naive],
    [`polycom_iterative`][poly.functional.polycom_iterative],
    [`polycom_horner`][poly.functional.polycom_horner]
    - for $q=x-s$: [`polyshift`][poly.functional.polyshift]
    - based on: [`polyval`][poly.functional.polyval]
    """
    q = tuple(q)
    match method:
        case 'naive':
            return polycom_naive(p, q)
        case 'iterative':
            return polycom_iterative(p, q)
        case 'horner':
            return polycom_horner(p, q)
        case _:
            raise ValueError('Invalid method')

def polyshift(p, s):
    """Return the polynomial `p` shifted by `s` on the abscissa.
    
    $$
        p(x - s)
    $$
    
    TODO: https://math.stackexchange.com/a/694571/1170417
    
    See also
    --------
    - based on: [`polycom`][poly.functional.polycom]
    """
    return polycom(p, polyaddc(polyx, -s))


#arithmetic
def polypos(p):
    """Return the polynomial with the unary positive operator applied.
    
    $$
        +p
    $$
    
    See also
    --------
    - wraps: [`vector.vecpos`](https://goessl.github.io/vector/functional/#vector.functional.vecpos)
    """
    return vecpos(p)

def polyneg(p):
    """Return the polynomial with the unary negative operator applied.
    
    $$
        -p
    $$
    
    See also
    --------
    - wraps: [`vector.vecneg`](https://goessl.github.io/vector/functional/#vector.functional.vecneg)
    """
    return vecneg(p)

def polyadd(*ps):
    r"""Return the sum of polynomials.
    
    $$
        p_0 + p_1 + \cdots
    $$
    
    See also
    --------
    - for constant or monomial: [`polyaddc`][poly.functional.polyaddc]
    - wraps: [`vector.vecadd`](https://goessl.github.io/vector/functional/#vector.functional.vecadd)
    """
    return vecadd(*ps)

def polyaddc(p, c, n=0):
    """Return the sum of polynomial `p` and a monomial of degree `n`.
    
    $$
        p + cx^n
    $$
    
    More efficient than `polyadd(p, polymonom(n, c))`.
    
    See also
    --------
    - for polynomial summand: [`polyadd`][poly.functional.polyadd]
    - wraps: [`vector.vecaddc`](https://goessl.github.io/vector/functional/#vector.functional.vecaddc)
    """
    return vecaddc(p, c, n)

def polysub(p, q):
    """Return the difference of two polynomials.
    
    $$
        p - q
    $$
    
    See also
    --------
    - for constant or monomial: [`polyaddc`][poly.functional.polyaddc]
    - wraps: [`vector.vecsub`](https://goessl.github.io/vector/functional/#vector.functional.vecsub)
    """
    return vecsub(p, q)

def polyscalarmul(a, p):
    r"""Return the product of a scalar and a polynomial.
    
    $$
        a\cdot p
    $$
    
    See also
    --------
    - for polynomial factor: [`polymul`][poly.functional.polymul]
    - wraps: [`vector.vecmul`](https://goessl.github.io/vector/functional/#vector.functional.vecmul)
    """
    return vecmul(a, p)

def polyscalartruediv(p, a):
    r"""Return the true division of a polynomial and a scalar.
    
    $$
        \frac{p}{a}
    $$
    
    See also
    --------
    - wraps: [`vector.vectruediv`](https://goessl.github.io/vector/functional/#vector.functional.vectruediv)
    """
    return vectruediv(p, a)

def polyscalarfloordiv(p, a):
    r"""Return the floor division of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor x^k
    $$
    
    See also
    --------
    - included in: [`polyscalardivmod`][poly.functional.polyscalardivmod]
    - wraps: [`vector.vecfloordiv`](https://goessl.github.io/vector/functional/#vector.functional.vecfloordiv)
    """
    return vecfloordiv(p, a)

def polyscalarmod(p, a):
    r"""Return the elementwise mod of a polynomial and a scalar.
    
    $$
        \sum_k\left(a_k \mod a \right)x^k
    $$
    
    See also
    --------
    - included in: [`polyscalardivmod`][poly.functional.polyscalardivmod]
    - wraps: [`vector.vecmod`](https://goessl.github.io/vector/functional/#vector.functional.vecmod)
    """
    return vecmod(p, a)

def polyscalardivmod(p, a):
    r"""Return the elementwise divmod of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor x^k, \ \sum_k\left(a_k \mod a \right)x^k
    $$
    
    See also
    --------
    - combines: [`polyscalarfloordiv`][poly.functional.polyscalarfloordiv], [`polyscalarmod`][poly.functional.polyscalarmod]
    - wraps: [`vector.vecdivmod`](https://goessl.github.io/vector/functional/#vector.functional.vecdivmod)
    """
    return vecdivmod(p, a)

def polymul_naive(p, q):
    """Return the product of two polynomials.
    
    $$
        p q
    $$
    
    Uses naive multiplication and summation.
    
    Both arguments must be sequences.
    
    See also
    --------
    - for any implementation: [`polymul`][poly.functional.polymul]
    - other implementations:
    [`polymul_karatsuba`][poly.functional.polymul_karatsuba]
    """
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
    if len(p) == 1:
        return (p[0] * q[0],)
    
    m = len(p) // 2
    pl, pu = p[:m], p[m:]
    ql, qu = q[:m], q[m:]
    
    z0 = _polymul_karatsuba(pl, ql)
    z2 = _polymul_karatsuba(pu, qu)
    z1 = _polymul_karatsuba(polyadd(pl, pu), polyadd(ql, qu))
    mid = polysub(polysub(z1, z0), z2)
    
    #return vecadd(z0, (0,)*m+mid, (0,)*(2*m)+z2)
    return z0[:m] + polyadd(z0[m:2*m], mid[:m]) + polyadd(z0[2*m:], mid[m:], z2)

def polymul_karatsuba(p, q):
    """Return the product of two polynomials.
    
    $$
        p q
    $$
    
    Uses the Karatsuba algorithm.
    
    Both arguments must be sequences.
    
    See also
    --------
    - for any implementation: [`polymul`][poly.functional.polymul]
    - other implementations:
    [`polymul_naive`][poly.functional.polymul_naive]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Karatsuba_algorithm>
    """
    if not p or not q:
        return ()
    
    p, q = sorted((q, p), key=len) #p shorter, q longer
    ql, qu = q[:len(p)], q[len(p):]
    rl, ru = _polymul_karatsuba(p, ql), polymul_naive(p, qu)
    #return vecadd(rl, (0,)*len(p)+ru)
    return rl[:len(p)] + polyadd(rl[len(p):], ru)

def polymul(*ps, method='naive'):
    r"""Return the product of polynomials.
    
    $$
        p_0 p_1 \cdots
    $$
    
    Available methods are
    
    - [`naive`][poly.functional.polymul_naive] &
    - [`karatsuba`][poly.functional.polymul_karatsuba].
    
    See also
    --------
    - implementations: [`polymul_naive`][poly.functional.polymul_naive],
    [`polymul_karatsuba`][poly.functional.polymul_karatsuba]
    - for scalar factor: [`polyscalarmul`][poly.functional.polyscalarmul]
    - for monomial factor: [`polymulx`][poly.functional.polymulx]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Karatsuba_algorithm>
    """
    ps = tuple(map(tuple, ps))
    match method:
        case 'naive':
            return reduce(polymul_naive, ps, polyone)
        case 'karatsuba':
            return reduce(polymul_karatsuba, ps, polyone)
        case _:
            raise ValueError('Invalid method')

def polymulx(p, n=1):
    """Return the product of polynomial `p` and a monomial of degree `n`.
    
    $$
        px^n
    $$
    
    More efficient than `polymul(p, polymonom(n, c))`.
    
    See also
    --------
    - for polynomial factor: [`polymul`][poly.functional.polymul]
    - wraps: [`vector.vecrshift`](https://goessl.github.io/vector/functional/#vector.functional.vecrshift)
    """
    return vecrshift(p, n)

def polypow_naive(p, n):
    """Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Uses repeated [multiplication][poly.functional.polymul].
    
    `p` must be a sequence.
    
    See also
    --------
    - for any implementation: [`polypow`][poly.functional.polypow]
    - other implementations:
    [`polypow_binary`][poly.functional.polypow_binary]
    """
    return reduce(polymul, repeat(p, n), polyone)

def polypow_binary(p, n):
    """Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Uses exponentiation by squaring.
    
    `p` must be a sequence.
    
    See also
    --------
    - for any implementation: [`polypow`][poly.functional.polypow]
    - other implementations:
    [`polypow_naive`][poly.functional.polypow_naive]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Exponentiation_by_squaring>
    """
    r = polyone
    while n:
        if n % 2 == 1:
            r = polymul(r, p)
        p = polymul(p, p)
        n //= 2
    return r

def polypow(p, n, method='binary'):
    """Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Available methods are 
    
    - [`naive`][poly.functional.polypow_naive] &
    - [`binary`][poly.functional.polypow_binary].
    
    See also
    --------
    - implementatiosn: [`polypow_naive`][poly.functional.polypow_naive],
    [`polypow_binary`][poly.functional.polypow_binary]
    
    References
    ----------
    - <https://en.wikipedia.org/wiki/Exponentiation_by_squaring>
    """
    p = tuple(p)
    match method:
        case 'naive':
            return polypow_naive(p, n)
        case 'binary':
            return polypow_binary(p, n)
        case _:
            raise ValueError('Invalid method')

#def polydiv(n, d):
#    """Return $(q, r)$ such that $n=qd+r$."""
#    #https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
#    n, d, q = vectrim(n), vectrim(d), polyzero
#    while n and polydeg(n)>=polydeg(d):
#        t = vecbasis(len(n)-len(d), n[-1]/d[-1])
#        q = vecadd(q, t)
#        n = vectrim(vecsub(n, polymul(d, t))[:-1])
#    return q, n


#calculus
def polyder(p, n=1):
    """Return the `n`-th derivative of polynomial `p`.
    
    $$
        p^{(n)}
    $$
    
    See also
    --------
    - inverse: [`polyantider`][poly.functional.polyantider]
    """
    #https://docs.python.org/3/library/itertools.html#itertools-recipes
    for _ in range(n):
        p = vechadamard(count(1), islice(p, 1, None))
    return p

def polyantider(p, n=1, c=0):
    """Return the `n`-th antiderivative of polynomial `p`.
    
    $$
        p^{(-n)}
    $$
    
    The integration constant `c` may be a sequence,
    or a scalar that is used in every integration step.
    
    Notes
    -----
    Integration is called antiderivative (`antider`)
    instead of integrate ('int') to avoid keyword collisions.
    For example `Poly((1, 2, 3)).int()`.
    
    See also
    --------
    - inverse: [`polyder`][poly.functional.polyder]
    """
    try:
        for i in range(n):
            p = (c[i],) + vechadamardtruediv(p, count(1))
        return p
    except TypeError:
        return polyantider(p, n, (c,)*n)


#sympy
def polysympify(p, symbol=spx):
    """Return the coefficient iterable `p` as a [`sympy.Poly`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly).
    
    See also
    --------
    - inverse: [`polyunsympify`][poly.functional.polyunsympify]
    """
    return sp.Poly.from_list(tuple(reversed(tuple(p))), symbol)

def polyunsympify(p):
    """Return [`sympy.Poly(p)`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly) as a coefficient `tuple`.
    
    See also
    --------
    - inverse: [`polysympify`][poly.functional.polysympify]
    """
    return tuple(reversed(p.all_coeffs()))
