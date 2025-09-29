from math import prod, factorial, comb
from operator import neg, mul
from fractions import Fraction
from functools import reduce, cache
from itertools import repeat, chain, islice, count, starmap, tee, zip_longest
from vector import *
import sympy as sp
from sympy.abc import x as spx

from .functional import *
from vector import *



__all__ = (#creation
           'H0', 'H1', 'H2',
           'herm', 'herm_recursive', 'herm_iterative', 'herm_explicit', 'hermgen',
           'hermzero', 'hermone', 'hermx',
           'hermmono', 'hermmono_recursive', 'hermmono_iterative', 'hermmono_explicit', 'hermmonogen',
           'hermrand', 'hermrandn', 'hermfromroots',
           #utility
           'hermeq', 'hermtrim', 'hermround', 'hermdeg',
           #conversion
           'herm2poly', 'herm2poly_naive', 'herm2poly_iterative', 'herm2poly_clenshaw',
           'poly2herm', 'poly2herm_naive', 'poly2herm_iterative', 'poly2herm_horner',
           #evaluation
           'hermval', 'hermval_naive', 'hermval_iterative', 'hermval_clenshaw',
           'hermvalzerogen', 'hermvalzero',
           #arithmetic
           'hermpos', 'hermneg', 'hermadd', 'hermaddc', 'hermsub',
           'hermscalarmul', 'hermscalartruediv', 'hermscalarfloordiv',
           'hermscalarmod', 'hermscalardivmod',
           'hermmul', 'hermmul_naive', 'hermmulx',
           'hermpow', 'hermpow_naive', 'hermpow_binary', 'hermmulpow',
           #calculus
           'hermder', 'hermantider',
           #sympy
           'hermsympify')



#creation
H0 = (1,)
r"""Zero-th Hermite polynomial in standard monomial basis.

$$
    H_0^{(P)}=\begin{pmatrix} 1 \end{pmatrix} \qquad H_0(x)=1
$$

A tuple with a single: `(1,)`.
"""
H1 = (0, 2)
r"""First Hermite polynomial in standard monomial basis.

$$
    H_0^{(P)}=\begin{pmatrix} 0 \\ 2 \end{pmatrix} \qquad H_1(x)=2x
$$

A tuple with a zero and a two: `(0, 2)`.
"""
H2 = (-2, 0, 4)
r"""Second Hermite polynomial in standard monomial basis.

$$
    H_2^{(P)}=\begin{pmatrix} -2 \\ 0 \\ 4 \end{pmatrix} \qquad H_2^{(P)}=4x^2-2
$$

A tuple with a negative two, a zero and a four: `(-2, 0, 4)`.
"""

def herm(n, method='iterative'):
    r"""Return the `n`-th Hermite polynomial in monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Available methods are
    
    - [`recursive`][poly.hermite_functional.herm_recursive],
    - [`iterative`][poly.hermite_functional.herm_iterative] &
    - [`explicit`][poly.hermite_functional.herm_explicit].
    
    Notes
    -----
    [`numpy.polynomial.hermite.herm2poly(vecbasis(n))`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.herm2poly.html) differs at $n\geq29$.
    
    See also
    --------
    - implementations: [`herm_recursive`][poly.hermite_functional.herm_recursive],
    [`herm_iterative`][poly.hermite_functional.herm_iterative],
    [`herm_explicit`][poly.hermite_functional.herm_explicit]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Explicit expression](https://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression)
    """
    match method:
        case 'recursive':
            return herm_recursive(n)
        case 'iterative':
            return herm_iterative(n)
        case 'explicit':
            return herm_explicit(n)
        case _:
            raise ValueError('Invalid method')

@cache
def herm_recursive(n):
    """Return the `n`-th Hermite polynomial in monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ recursively.
    
    Cached, therefore fast if used repeatedly.
    
    See also
    --------
    - for any implementation: [`herm`][poly.hermite_functional.herm]
    - other implementations: [`herm_iterative`][poly.hermite_functional.herm_iterative],
    [`herm_explicit`][poly.hermite_functional.herm_explicit]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    if n == 0:
        return H0
    elif n == 1:
        return H1
    else:
        return polyscalarmul(2, polysub(polymulx(herm_recursive(n-1)),
                                        polyscalarmul(n-1, herm_recursive(n-2))))

def herm_iterative(n):
    """Return the `n`-th Hermite polynomial in monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ iteratively.
    
    See also
    --------
    - for any implementation: [`herm`][poly.hermite_functional.herm]
    - other implementations: [`herm_recursive`][poly.hermite_functional.herm_recursive],
    [`herm_explicit`][poly.hermite_functional.herm_explicit]
    - uses: [`hermgen`][poly.hermite_functional.hermgen]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    return next(islice(hermgen(), n, None))

def herm_explicit(n):
    r"""Return the `n`-th Hermite polynomial in monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Using the explicit expression $H_n(x)=n!\sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}\frac{(-1)^m}{m!(n-2m)!}(2x)^{n-2m}$.
    
    Notes
    -----
    The leading coefficient of the $n$-th Hermite polynomial is $2^n$.
    
    Every second coefficient is zero.
    
    Every other one is set by $n!\sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}\frac{(-1)^m}{m!(n-2m)!}(2x)^{n-2m}$ resulting in
    
    $$
        \begin{aligned}
            a_{n, n-2m} &= n!2^{n-2m}\frac{(-1)^m}{m!(n-2m)!} &&\mid i=n-2m \Leftrightarrow m=\frac{n-i}{2} \\
            a_{n, i} &= n!2^i\frac{(-1)^\frac{n-i}{2}}{\frac{n-i}{2}!i!}
        \end{aligned}
    $$
    
    and following the recursion relation
    
    $$
        \begin{aligned}
            \frac{a_{n, i}}{a_{n, i+2}} &= \frac{n!2^i\frac{(-1)^\frac{n-i}{2}}{\frac{n-i}{2}!i!}}{n!2^{i+2}\frac{(-1)^\frac{n-i-2}{2}}{\frac{n-i-2}{2}!(i+2)!}} \\
            &= \frac{n!2^i(-1)^\frac{n-i}{2}\frac{n-i-2}{2}!(i+2)!}{n!2^{i+2}(-1)^\frac{n-i-2}{2}\frac{n-i}{2}!i!} \\
            &= -\frac{(i+1)(i+2)}{2(n-i)} \\
            &= -\frac{(n-2m+1)(n-2m+2)}{2(n-n+2m)} \\
            &= -\frac{(n-2m+1)(n-2m+2)}{4m}
        \end{aligned}
    $$
    
    which allows setting the coefficients in descending order iteratively without factorials.
    
    See also
    --------
    - for any implementation: [`herm`][poly.hermite_functional.herm]
    - other implementations: [`herm_recursive`][poly.hermite_functional.herm_recursive],
    [`herm_iterative`][poly.hermite_functional.herm_iterative]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Explicit expression](https://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression)
    """
    r = [0] * (n+1)
    r[-1] = 2**n
    for m in range(1, n//2+1):
        r[n-2*m] = -r[n-2*(m-1)] * (n-2*m+1)*(n-2*m+2) // (4*m)
    return tuple(r)

def hermgen():
    r"""Yield the Hermite polynomials in standard monomial basis.
    
    $$
    (H_0^{(P)})_{n\in\mathbb{N}_0} = (H_0^{(P)}, H_2^{(P)}, H_2^{(P)}, \dots)
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ iteratively.
    
    See also
    --------
    - used in: [`herm_iterative`][poly.hermite_functional.herm_iterative]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    yield (Hkm1 := H0)
    yield (Hk := H1)
    for km1 in count(1):
        Hkm1, Hk = Hk, polyscalarmul(2, polysub(polymulx(Hk), polyscalarmul(km1, Hkm1)))
        yield Hk

hermzero = veczero
"""Zero Hermite polynomial series.

$$
    0^{(H)}
$$

An empty tuple: `()`.

See also
--------
- in standard monomial basis: [`polyzero`][poly.functional.polyzero]
- wraps: [`vector.veczero`](https://goessl.github.io/vector/functional/#vector.functional.veczero)

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermzero`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermzero.html)
"""
hermone = (1, )
r"""Constant one Hermite polynomial series.

$$
    1^{(H)}=\begin{pmatrix}1\end{pmatrix} \qquad 1=H_0
$$

A tuple with a single one: `(1,)`.

Notes
-----
A `Fraction` `tuple` (`(Fraction(1), )`) would be more consitent with [`hermx`][poly.hermite_functional.hermx], but would then require `Fraction` arithmetic for multiplicative functions, ([`hermmul`][poly.hermite_functional.hermmul], [`hermpow`][poly.hermite_functional.hermpow]).

See also
--------
- in standard monomial basis: [`polyone`][poly.functional.polyone]

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermone`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermone.html)
"""
hermx = (Fraction(0), Fraction(1, 2))
r"""Identity Hermite polynomial series.

$$
    x^{(H)}=\begin{pmatrix} 0 \\ \frac{1}{2} \end{pmatrix} \qquad x=\frac{1}{2}H_1(x)
$$

A tuple with a zero and a half: `(0, 1/2)`.

See also
--------
- in standard monomial basis: [`polyx`][poly.functional.polyx]

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermx.html)
"""

def hermmono(n, method='explicit'):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        (x^n)^{(H)}
    $$
    
    Available methods are:
    
    - [`recursive`][poly.hermite_functional.hermmono_recursive],
    - [`iterative`][poly.hermite_functional.hermmono_iterative] &
    - [`explicit`][poly.hermite_functional.hermmono_explicit].
    
    See also
    --------
    - implementations: [`hermmono_recursive`][poly.hermite_functional.hermmono_recursive],
    [`hermmono_iterative`][poly.hermite_functional.hermmono_iterative],
    [`hermmono_explicit`][poly.hermite_functional.hermmono_explicit]
    - for all $n$: [`hermmonogen`][poly.hermite_functional.hermmonogen]
    """
    match method:
        case 'recursive':
            return hermmono_recursive(n)
        case 'iterative':
            return hermmono_iterative(n)
        case 'explicit':
            return hermmono_explicit(n)
        case _:
            raise ValueError('Invalid method')

@cache
def hermmono_recursive(n):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        (x^n)^{(H)}
    $$
    
    Using [`hermmulx`][poly.hermite_functional.hermmulx] recursively.
    
    Cached, therefore fast if used repeatedly.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite_functional.hermmono]
    - other implementations: [`hermmono_iterative`][poly.hermite_functional.hermmono_iterative],
    [`hermmono_explicit`][poly.hermite_functional.hermmono_explicit]
    """
    if n == 0:
        return (Fraction(1),) #hermone but as `Fraction`
    return hermmulx(hermmono_recursive(n-1))

@cache
def hermmono_iterative(n):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        (x^n)^{(H)}
    $$
    
    Using [`hermmulx`][poly.hermite_functional.hermmulx] iteratively.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite_functional.hermmono]
    - other implementations: [`hermmono_recursive`][poly.hermite_functional.hermmono_recursive],
    [`hermmono_explicit`][poly.hermite_functional.hermmono_explicit]
    - uses: [`hermmonogen`][poly.hermite_functional.hermmonogen]
    """
    return next(islice(hermmonogen(), n, None))

def hermmono_explicit(n):
    r"""Return `x^n` as Hermite polynomial series.
    
    $$
        (x^n)^{(H)}
    $$
        
    Using the explicit expression $x^n=\frac{n!}{2^n}\sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}\frac{1}{m!(n-2m)!}H_{n-2m}(x)$.
    
    Notes
    -----
    The leading coefficient is $\frac{1}{2^n}$.
    
    Every second coefficient is zero.
    
    Every other one is set by $x^n=\frac{n!}{2^n}\sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}\frac{H_{n-2m}(x)}{m!(n-2m)!}$ resulting in
    
    $$
        \begin{aligned}
            a^{-1}_{n, n-2m} &= \frac{n!}{2^nm!(n-2m)!} &&\mid i=n-2m \Leftrightarrow m=\frac{n-i}{2} \\
            a^{-1}_{n, i} &= \frac{n!}{2^n\frac{n-i}{2}!i!}
        \end{aligned}
    $$
    
    and following the recursion relation
    
    $$
        \begin{aligned}
            \frac{a_{n, i}}{a_{n, i+2}} &= \frac{\frac{n!}{2^n\frac{n-i}{2}!i!}}{\frac{n!}{2^n\frac{n-i-2}{2}!(i+2)!}} \\
            &= \frac{n!2^n\frac{n-i-2}{2}!(i+2)!}{n!2^n\frac{n-i}{2}!i!} \\
            &= \frac{2(i+1)(i+2))}{n-i} \\
            &= \frac{2(n-2m+1)(n-2m+2))}{n-n-2m} \\
            &= -\frac{(n-2m+1)(n-2m+2))}{m}
        \end{aligned}
    $$
    
    which allows setting the coefficients in descending order iteratively without factorials.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite_functional.hermmono]
    - other implementations: [`hermmono_recursive`][poly.hermite_functional.hermmono_recursive],
    [`hermmono_explicit`][poly.hermite_functional.hermmono_explicit]
    - uses: [`hermmonogen`][poly.hermite_functional.hermmonogen]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Inverse explicit expression](https://en.wikipedia.org/wiki/Hermite_polynomials#Inverse_explicit_expression)
    """
    r = [Fraction()] * (n + 1)
    r[n] = Fraction(1, 2**n)
    for m in range(1, n//2+1):
        r[n-2*m] = r[n-2*(m-1)] * ((n-2*m+1)*(n-2*m+2)) / m
    return tuple(r)

def hermmonogen():
    r"""Yield standard monomials `x^n` as Hermite polynomial series.
    
    $$
        \left((x^n)^{(H)}\right)_{n\in\mathbb{N}_0} = (1^{(H)}, x^{(H)}, (x^2)^{(H)}, \dots)
    $$
    
    Using [`hermmulx`][poly.hermite_functional.hermmulx] iteratively.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - used in: [`hermmono_iterative`][poly.hermite_functional.hermmono_iterative]
    """
    yield (h := (Fraction(1),)) #hermone but as `Fraction`
    for i in count():
        h = hermmulx(h)
        yield h

def hermrand(n):
    r"""Return a random Hermite polynomial series of degree `n`.
    
    $$
        \sum_{k=0}^na_kH_k(x) \quad a_k\sim\mathcal{U}([0, 1[)
    $$
    
    The coefficients are sampled from a uniform distribution in `[0, 1[`.
    
    See also
    --------
    - in standard monomial basis: [`polyrand`][poly.functional.polyrand]
    - wraps: [`vector.vecrand`](https://goessl.github.io/vector/functional/#vector.functional.vecrand)
    """
    return vecrand(n+1)

def hermrandn(n, normed=True, mu=0, sigma=1):
    r"""Return a random Hermite polynomial series of degree `n`.
    
    $$
        \sum_{k=0}^na_kH_k(x) \quad a_k\sim\mathcal{N}(\mu, \sigma)
    $$
    
    The coefficients are sampled from a normal distribution.
    
    Normed with respect to the vector norm $\sum_k|a_k|^2$.
    TODO: normalise with respect to Hermite polynomial weighted norm.
    
    See also
    --------
    - in standard monomial basis: [`polyrandn`][poly.functional.polyrandn]
    - wraps: [`vector.vecrandn`](https://goessl.github.io/vector/functional/#vector.functional.vecrandn)
    """
    return vecrandn(n+1, normed=normed, mu=mu, sigma=sigma)

def hermfromroots(*xs):
    r"""Return the Hermite polynomial series with the given roots.
    
    $$
        \prod_k\left(\frac{H_1}{2}-x_k\right)
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyfromroots`][poly.functional.polyfromroots]
    """
    return hermmul(*zip(map(neg, xs), repeat(Fraction(1, 2))))


#utility
def hermeq(g, h):
    r"""Return if two Hermite polynomial series are equal.
    
    $$
        g \overset{?}{=} h
    $$
    
    See also
    --------
    - standard monomial basis: [`polyeq`][poly.functional.polyeq]
    - wraps: [`vector.veceq`](https://goessl.github.io/vector/functional/#vector.functional.veceq)
    """
    return veceq(g, h)

def hermtrim(h, tol=1e-9):
    """Remove all leading near zero (`abs(h_i)<=tol`) coefficients.
    
    See also
    --------
    - standard monomial basis: [`polytrim`][poly.functional.polytrim]
    - wraps: [`vector.vectrim`](https://goessl.github.io/vector/functional/#vector.functional.vectrim)
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermtrim`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermtrim.html)
    """
    return vectrim(h, tol)

def hermround(h, ndigits=None):
    r"""Round all coefficients to the given precision.
    
    $$
        \sum_k\text{round}_\text{ndigits}(h_k)H_k
    $$
    
    See also
    --------
    - standard monomial basis: [`polyround`][poly.functional.polyround]
    - wraps: [`vector.vecround`](https://goessl.github.io/vector/functional/#vector.functional.vecround)
    """
    return vecround(h, ndigits)

def hermdeg(h):
    r"""Return the degree of a Hermite polynomial series.
    
    $$
        \deg(h)
    $$
    
    $\deg(0)=-1$ is used for the empty zero Hermite polynomial series.
    
    Doesn't handle leading zeros, use
    [`hermtrim`][poly.hermite_functional.hermtrim] if needed.
    
    `h` must support `len(h)`.
    
    Notes
    -----
    $\deg(0)=-\infty$ is more commonly used but the expected return type is an
    `int` and `-math.inf` is of type `float`. Therefore the $\deg(0)=-1$
    convention was choosen to conserve the return type.
    
    See also
    --------
    - standard monomial basis: [`polydeg`][poly.functional.polydeg]
    
    References
    ----------
    - [Wikipedia - Degree of a polynomial - Degree of the zero polynomial](https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial)
    """
    return len(h) - 1


#conversion
def herm2poly(h, method='iterative'):
    """Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Available methods are 
    
    - [`naive`][poly.hermite_functional.herm2poly_naive],
    - [`iterative`][poly.hermite_functional.herm2poly_iterative] &
    - [`clenshaw`][poly.hermite_functional.herm2poly_clenshaw].
    
    See also
    --------
    - implementations: [`herm2poly_naive`][poly.hermite_functional.herm2poly_naive],
    [`herm2poly_iterative`][poly.hermite_functional.herm2poly_iterative],
    [`herm2poly_clenshaw`][poly.hermite_functional.herm2poly_clenshaw]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.herm2poly`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.herm2poly.html)
    """
    match method:
        case 'naive':
            return herm2poly_naive(h)
        case 'iterative':
            return herm2poly_iterative(h)
        case 'clenshaw':
            return herm2poly_clenshaw(h)
        case _:
            raise ValueError('Invalid method')

def herm2poly_naive(h):
    """Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Uses naive Hermite polynomial creation by the default method of
    [`herm`][poly.hermite_functional.herm].
    
    See also
    --------
    - for any implementation: [`herm2poly`][poly.hermite_functional.herm2poly]
    - other implementations: [`herm2poly_iterative`][poly.hermite_functional.herm2poly_iterative],
    [`herm2poly_clenshaw`][poly.hermite_functional.herm2poly_clenshaw]
    """
    return polyadd(*(polyscalarmul(hi, herm(i)) for i, hi in enumerate(h)))

def herm2poly_iterative(h):
    """Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Uses iterative Hermite polynomial creation like
    [`herm_iterative`][poly.hermite_functional.herm_iterative].
    
    See also
    --------
    - for any implementation: [`herm2poly`][poly.hermite_functional.herm2poly]
    - other implementations: [`herm2poly_naive`][poly.hermite_functional.herm2poly_naive],
    [`herm2poly_clenshaw`][poly.hermite_functional.herm2poly_clenshaw]
    - uses: [`hermgen`][poly.hermite_functional.hermgen]
    """
    return polyadd(*starmap(polyscalarmul, zip(h, hermgen())))

def herm2poly_clenshaw(h):
    """Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Uses the Clenshaw algorithm like
    [`hermval_clenshaw`][poly.hermite_functional.hermval_clenshaw].
    
    See also
    --------
    - for any implementation: [`herm2poly`][poly.hermite_functional.herm2poly]
    - other implementations [`herm2poly_naive`][poly.hermite_functional.herm2poly_naive],
    [`herm2poly_iterative`][poly.hermite_functional.herm2poly_iterative]
    """
    h = tuple(h)
    if not h:
        return hermzero
    a, b = hermzero, hermzero
    for n, hn in reversed(tuple(enumerate(h[1:], 1))):
        a, b = polyaddc(polysub(polyscalarmul(2, polymulx(a)), polyscalarmul(2*(n+1), b)), hn), a
    return polyaddc(polysub(polyscalarmul(2, polymulx(a)), polyscalarmul(2, b)), h[0])

def poly2herm(p, method='iterative'):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Available methods are 
    
    - [`naive`][poly.hermite_functional.poly2herm_naive],
    - [`iterative`][poly.hermite_functional.poly2herm_iterative] &
    - [`horner`][poly.hermite_functional.poly2herm_horner].
    
    See also
    --------
    - implementations: [`poly2herm_naive`][poly.hermite_functional.poly2herm_naive],
    [`poly2herm_iterative`][poly.hermite_functional.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite_functional.poly2herm_horner]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.poly2herm`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.poly2herm.html)
    """
    match method:
        case 'naive':
            return poly2herm_naive(p)
        case 'iterative':
            return poly2herm_iterative(p)
        case 'horner':
            return poly2herm_horner(p)
        case _:
            raise ValueError('Invalid method')

def poly2herm_naive(p):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Uses naive monomial as Hermite polynomial series creation by the default
    method of [`hermmono`][poly.hermite_functional.hermmono].
    
    See also
    --------
    - any implementation: [`poly2herm`][poly.hermite_functional.poly2herm]
    - other implementations: [`poly2herm_iterative`][poly.hermite_functional.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite_functional.poly2herm_horner]
    """
    return hermadd(*(hermscalarmul(pi, hermmono(i)) for i, pi in enumerate(p)))

def poly2herm_iterative(p):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Uses iterative monomial creation of [`hermmonogen`][poly.hermite_functional.hermmonogen].
    
    See also
    --------
    - for any implementation: [`poly2herm`][poly.hermite_functional.poly2herm]
    - other implementations: [`poly2herm_iterative`][poly.hermite_functional.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite_functional.poly2herm_horner]
    """
    return hermadd(*starmap(hermscalarmul, zip(p, hermmonogen())))

def poly2herm_horner(p):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Uses Horner's method.
    
    `p` must be reversible.
    
    See also
    --------
    - for any implementation: [`poly2herm`][poly.hermite_functional.poly2herm]
    - other implementations: [`poly2herm_iterative`][poly.hermite_functional.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite_functional.poly2herm_horner]
    
    References
    ----------
    - [Wikipedia - Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method)
    """
    r = hermzero
    for pi in reversed(p):
        r = hermaddc(hermmulx(r), pi)
    return r


#evaluation
def hermval(h, x, method='clenshaw'):
    """Return the value of Hermite polynomial series `h` evaluated at point `x`.
    
    $$
        h(x)
    $$
    
    Available methods are
    
    - [`naive`][poly.hermite_functional.hermval_naive],
    - [`iterative`][poly.hermite_functional.hermval_iterative] &
    - [`clenshaw`][poly.hermite_functional.hermval_clenshaw].
    
    See also
    --------
    - implementations: [`hermval_naive`][poly.hermite_functional.hermval_naive],
    [`hermval_iterative`][poly.hermite_functional.hermval_iterative],
    [`hermval_clenshaw`][poly.hermite_functional.hermval_clenshaw]
    - in standard monomial basis: [`polyval`][poly.functional.polyval]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermval`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermval.html)
    """
    match method:
        case 'naive':
            return hermval_naive(h, x)
        case 'iterative':
            return hermval_iterative(h, x)
        case 'clenshaw':
            return hermval_clenshaw(h, x)
        case _:
            raise ValueError('Invalid method')

def hermval_naive(h, x):
    """Return the value of Hermite polynomial series `h` evaluated at point `x`.
    
    $$
        h(x)
    $$
    
    Converts $h$ to standard monomial basis and evaluates with [`polyval`][poly.functional.polyval].
    
    See also
    --------
    - for any implementation: [`hermval`][poly.hermite_functional.hermval]
    - other implementations: [`hermval_iterative`][poly.hermite_functional.hermval_iterative],
    [`hermval_clenshaw`][poly.hermite_functional.hermval_clenshaw]
    - in standard monomial basis: [`polyval_naive`][poly.functional.polyval_naive]
    """
    return polyval(herm2poly(h), x)

def hermval_iterative(h, x):
    """Return the value of Hermite polynomial series `h` evaluated at point `x`.
    
    $$
        h(x)
    $$
    
    Uses iterative `H_n(x)` generation.
    
    See also
    --------
    - for any implementation: [`hermval`][poly.hermite_functional.hermval]
    - other implementations: [`hermval_naive`][poly.hermite_functional.hermval_naive],
    [`hermval_clenshaw`][poly.hermite_functional.hermval_clenshaw]
    - uses: [`hermvalgen`][poly.hermite_functional.hermvalgen]
    - in standard monomial basis: [`polyval_iterative`][poly.functional.polyval_iterative]
    """
    return sum(map(mul, h, hermvalgen(x)), start=type(x)(0))

def hermval_clenshaw(h, x):
    """Return the value of Hermite polynomial series `h` evaluated at point `x`.
    
    $$
        h(x)
    $$
    
    Uses the Clenshaw algorithm.
    
    See also
    --------
    - for any implementation: [`hermval`][poly.hermite_functional.hermval]
    - other implementations: [`hermval_naive`][poly.hermite_functional.hermval_naive],
    [`hermval_iterative`][poly.hermite_functional.hermval_iterative]
    - in standard monomial basis: [`polyval_horner`][poly.functional.polyval_horner]
    
    References
    ----------
    - [Wikipedia - Clenshaw algorithm](https://en.wikipedia.org/wiki/Clenshaw_algorithm)
    """
    h = tuple(h)
    if not h:
        return type(x)(0)
    a, b = 0, 0
    for n, hn in reversed(tuple(enumerate(h[1:], 1))):
        a, b = hn + 2*x*a - 2*(n+1)*b, a
    return h[0] + 2*x*a - 2*b

def hermvalgen(x):
    r"""Yield the value of Hermite polynomials evaluated at point `x`.
    
    $$
        (H_n(x))_{n\in\mathbb{N}_0} = (H_0(x), H_1(x), H_2(x), \dots)
    $$
    
    Using the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ iteratively.
    
    See also
    --------
    - used in: [`hermval_iterative`][poly.hermite_functional.hermval_iterative]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    yield (Hkm1 := type(x)(1))
    yield (Hk := 2*x)
    for k in count(1):
        Hkm1, Hk = Hk, 2*(Hk*x - k*Hkm1)
        yield Hk

def hermvalzerogen():
    r"""Yield the values of Hermite polynomials evaluated at point 0.
    
    $$
        (H_n(0))_{n\in\mathbb{N}_0} = (H_0(0), H_1(0), H_2(0), \dots)
    $$
    
    Notes
    -----
    $$
        \begin{gathered}
            H_n(0) = \begin{cases}
                0 & n\in\mathbb{U} \\
                (-2)^\frac{n}{2}(n-1)!! & n\in\mathbb{G}
            \end{cases} \\
            1, 0, -2, 0, 12, 0, -120, 0, 1680, 0, -30240, \dots
        \end{gathered}
    $$
    
    See also
    --------
    - for any series: [`hermvalzero`][poly.hermite_functional.hermvalzero]
    """
    yield (y := 1)
    for i in count(1, 2):
        yield 0
        y *= -2 * i
        yield y

def hermvalzero(h):
    """Return the value of Hermite polynomial series `h` evaluated at point 0.
    
    $$
        h(0)
    $$
    
    More efficient than `hermval(h, 0)`.
    
    See also
    --------
    - for any point: [`hermval`][poly.hermite_functional.hermval]
    """
    return sum(map(prod, islice(zip(h, hermvalzerogen()), 0, None, 2)))


#arithmetic
def hermpos(h):
    """Return the Hermite polynomial series with the unary positive operator applied.
    
    $$
        +h
    $$
    
    See also
    --------
    - in standard monomial basis: [`polypos`][poly.functional.polypos]
    - wraps: [`vector.vecpos`](https://goessl.github.io/vector/functional/#vector.functional.vecpos)
    """
    return vecpos(h)

def hermneg(h):
    """Return the Hermite polynomial series with the unary negative operator applied.
    
    $$
        -h
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyneg`][poly.functional.polyneg]
    - wraps: [`vector.vecneg`](https://goessl.github.io/vector/functional/#vector.functional.vecneg)
    """
    return vecneg(h)

def hermadd(*hs):
    r"""Return the sum of Hermite polynomial series.
    
    $$
        h_0 + h_1 + \cdots
    $$
    
    See also
    --------
    - for single coefficient summand ($cH_n$): [`hermaddc`][poly.hermite_functional.hermaddc]
    - in standard monomial basis: [`polyadd`][poly.functional.polyadd]
    - wraps: [`vector.vecadd`](https://goessl.github.io/vector/functional/#vector.functional.vecadd)
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermadd`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermadd.html)
    """
    return vecadd(*hs)

def hermaddc(h, c, n=0):
    """Return the sum of Hermite polynomial series `h` and the `n`-th Hermite polynomial.
    
    $$
        h + cH_n
    $$
    
    More efficient than `hermadd(h, vecbasis(n, c))`.
    
    See also
    --------
    - for series summand: [`hermadd`][poly.hermite_functional.hermadd]
    - in standard monomial basis: [`polyaddc`][poly.functional.polyaddc]
    - wraps: [`vector.vecaddc`](https://goessl.github.io/vector/functional/#vector.functional.vecaddc)
    """
    return vecaddc(h, c, n)

def hermsub(g, h):
    """Return the difference of two Hermite polynomial series.
    
    $$
        g - h
    $$
    
    See also
    --------
    - for single coefficient subtrahend ($cH_n$): [`hermaddc`][poly.hermite_functional.hermaddc]
    - in standard monomial basis: [`polysub`][poly.functional.polysub]
    - wraps: [`vector.vecsub`](https://goessl.github.io/vector/functional/#vector.functional.vecsub)
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermsub`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermsub.html)
    """
    return vecsub(g, h)

def hermscalarmul(a, h):
    r"""Return the product of a scalar and a Hermite polynomial series.
    
    $$
        a\cdot h
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyscalarmul`][poly.functional.polyscalarmul]
    - wraps: [`vector.vecmul`](https://goessl.github.io/vector/functional/#vector.functional.vecmul)
    """
    return vecmul(a, h)

def hermscalartruediv(h, a):
    r"""Return the true division of a Hermite polynomial series and a scalar.
    
    $$
        \frac{h}{a}
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyscalartruediv`][poly.functional.polyscalartruediv]
    - wraps: [`vector.vectruediv`](https://goessl.github.io/vector/functional/#vector.functional.vectruediv)
    """
    return vectruediv(h, a)

def hermscalarfloordiv(h, a):
    r"""Return the floor division of a Hermite polynomial series and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor H_k
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyscalarfloordiv`][poly.functional.polyscalarfloordiv]
    - wraps: [`vector.vecmfloordiv`](https://goessl.github.io/vector/functional/#vector.functional.vecmfloordiv)
    """
    return vecfloordiv(h, a)

def hermscalarmod(h, a):
    r"""Return the elementwise mod of a Hermite polynomial series and a scalar.
    
    $$
        \sum_k\left(a_k \mod a \right)H_k
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyscalarmod`][poly.functional.polyscalarmod]
    - wraps: [`vector.vecmod`](https://goessl.github.io/vector/functional/#vector.functional.vecmod)
    """
    return vecmod(h, a)

def hermscalardivmod(h, a):
    r"""Return the elementwise divmod of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor H_k, \ \sum_k\left(a_k \mod a \right)H_k
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyscalardivmod`][poly.functional.polyscalardivmod]
    - wraps: [`vector.vecdivmod`](https://goessl.github.io/vector/functional/#vector.functional.vecdivmod)
    """
    return vecdivmod(h, a)

def hermmul(*hs, method='naive'):
    r"""Return the product of Hermite polynomial series.
    
    $$
        h_0 h_1 \cdots
    $$
    
    Available methods are
    
    - [`naive`][poly.hermite_functional.hermmul_naive].
    
    TODO: Karatsuba?
    
    See also
    --------
    - implementations: [`hermmul_naive`][poly.hermite_functional.hermmul_naive]
    - for monomial factor: [`hermmulx`][poly.hermite_functional.hermmulx]
    - in standard monomial basis: [`polymul`][poly.functional.polymul]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermmul`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermmul.html)
    """
    hs = tuple(map(tuple, hs))
    match method:
        case 'naive':
            return reduce(hermmul_naive, hs, hermone)
        case _:
            raise ValueError('Invalid method')

def hermmul_naive(g, h):
    r"""Return the product of two Hermite polynomials series.
    
    $$
        g h
    $$
    
    Uses naive multiplication and summation.
    
    Both arguments must be sequences.
    
    Notes
    -----
    From [1] eq. A.8 we know that
    
    $$
        \tilde{H}_i\tilde{H}_j = \sum_{k=0}^{\min\{i,j\}}k!\binom{i}{k}\binom{j}{k}\tilde{H}_{i+j-2k}.
    $$
    
    Therefore we have
    
    $$
        \begin{aligned}
            &H_i(x)H_j(x) &&\mid H_j(x) = 2^\frac{j}{2}\tilde{H}_j(\sqrt{2}x) \\
            &= 2^\frac{i+j}{2} \tilde{H}_i(\sqrt{2}x) \tilde{H}_j(\sqrt{2}x) &&\mid \tilde{H}_i\tilde{H}_j = \sum_{k=0}^{\min\{i, j\}}k!\binom{i}{k}\binom{j}{k}\tilde{H}_{i+j-2k} \\
            &= 2^\frac{i+j}{2} \sum_{k=0}^{\min\{i, j\}} k!\binom{i}{k}\binom{j}{k} \tilde{H}_{i+j-2k}(\sqrt{2}x) &&\mid \cdot1=2^{k-k} \\
            &= \sum_{k=0}^{\min\{i, j\}} 2^kk!\binom{i}{k}\binom{j}{k}2^\frac{i+j-2k}{2} \tilde{H}_{i+j-2k}(\sqrt{2}x) &&\mid H_j(x)=2^\frac{j}{2} \tilde{H}_j(\sqrt{2}x) \\
            &= \sum_{k=0}^{\min\{i, j\}}2^kk!\binom{i}{k}\binom{j}{k} H_{i+j-2k}(x)
        \end{aligned}
    $$
    
    See also
    --------
    - for any implementation: [`hermmul`][poly.hermite_functional.hermmul]
    - in standard monomial basis: [`polymul_naive`][poly.functional.polymul_naive]
    
    References
    ----------
    - [1] Zhi-yuan Huang & Jia-an Yan: Introduction to Infinite Dimensional Stochastic Analysis. DOI: 10.1007/978-94-011-4108-6
    - [Wikipedia - Hermite polynomials - Definition](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition)
    """
    g, h = tuple(g), tuple(h)
    if not g or not h:
        return hermzero
    r = [0] * (len(g)+len(h)-1)
    for i, gi in enumerate(g):
        for j, hj in enumerate(h):
            for k in range(min(i, j)+1):
                r[i+j-2*k] += 2**k*factorial(k)*comb(i, k)*comb(j, k) * gi*hj
    return tuple(r)

def hermmulx(h):
    r"""Return the product of Hermite polynomial series `h` and a monomial of degree `n`.
    
    $$
        h(x)x^n
    $$
    
    More efficient than `hermmul(h, hermx)`.
    
    Notes
    --------
    For a single Hermite polynomial
    
    $$
        \begin{aligned}
            H_{k+1}(x) &= 2xH_k(x)-H_k'(x) \\
            H_k'(x) &= 2kH_{k-1}(x) \\
            \Rightarrow \quad xH_k(x) &= \frac{H_{k+1}(x)}{2}+kH_{k-1}(x)
        \end{aligned}
    $$
    
    and therefore for a series
    
    $$
        \begin{aligned}
            xh(x) &= x\sum_kh_kH_k(x) \\
            &= \sum_kh_kxH_k(x) \\
            &= \sum_kh_k\left(\frac{1}{2}H_{k+1}(x)+kH_{k-1}(x)\right) \\
            &= \sum_k\frac{h_k}{2}H_{k+1}(x)+\sum_kkh_kH_{k-1}(x) \\
            &= \sum_k\frac{h_{k-1}}{2}H_k(x)+\sum_k(k+1)h_{k+1}H_k(x) \\
            &= \sum_k\left(\frac{h_{k-1}}{2}+(k+1)h_{k+1}\right)H_k(x)
        \end{aligned}
    $$
    
    See also
    --------
    - for Hermite polynomial series as factor: [`hermmul`][poly.hermite_functional.hermmul]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermmulx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermmulx.html)
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    im1, ip1 = tee(h)
    _, h1 = next(ip1, None), next(ip1, 0)
    return tuple(chain((h1,), (hkm1/2+(k+1)*hkp1 for k, (hkm1, hkp1) in
            enumerate(zip_longest(im1, ip1, fillvalue=0), 1))))

def hermpow(h, n, method='binary'):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    Available methods are 
    
    - [`naive`][poly.hermite_functional.hermpow_naive] &
    - [`binary`][poly.hermite_functional.hermpow_binary].
    
    See also
    --------
    - implementations: [`hermpow_naive`][poly.hermite_functional.hermpow_naive],
    [`hermpow_binary`][poly.hermite_functional.hermpow_binary],
    - in standard monomial basis: [`polypow`][poly.functional.polypow]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermpow`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermpow.html)
    """
    match method:
        case 'naive':
            return hermpow_naive(h, n)
        case 'binary':
            return hermpow_binary(h, n)
        case _:
            raise ValueError('Invalid method')

def hermpow_naive(h, n):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    Uses repeated multiplication.
    
    `h` must be a sequence.
    
    See also
    --------
    - for any implementation: [`hermpow`][poly.hermite_functional.hermpow]
    - other implementations: [`hermpow_binary`][poly.hermite_functional.hermpow_binary]
    - in standard monomial basis: [`polypow_naive`][poly.functional.polypow_naive]
    """
    return reduce(hermmul, repeat(h, n), hermone)

def hermpow_binary(h, n):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    Uses exponentiation by squaring.
    
    `h` must be a sequence.
    
    See also
    --------
    - for any implementations: [`hermpow`][poly.hermite_functional.hermpow]
    - other implementations[`hermpow_naive`][poly.hermite_functional.hermpow_naive]
    - in standard monomial basis: [`polypow_binary`][poly.functional.polypow_binary]
    
    References
    ----------
    - [Wikipedia - Exponentiation by squaring](https://en.wikipedia.org/wiki/Exponentiation_by_squaring)
    """
    r = hermone
    while n:
        if n % 2 == 1:
            r = hermmul(r, h)
        h = hermmul(h, h)
        n //= 2
    return r

def hermmulpow(alpha):
    r"""Return product of Hermite Polynomials.

    $$
        H^\alpha = \prod_iH_i^{\alpha_i}
    $$
    """
    r = hermone
    for i, ai in enumerate(alpha):
        r = hermmul(r, hermpow(vecbasis(i), ai))
    return r


#calculus
def hermder(h, n=1):
    r"""Return the `n`-th derivative of Hermite polynomial series `h`.
    
    $$
        h^{(n)}
    $$
    
    Notes
    -----
    $$
        \begin{aligned}
            H_k' &= 2kH_{k-1} \\
            h'(x) &= \frac{d}{dx}\sum_kh_kH_k(x) \\
            &= \sum_kh_kH_k'(x) \\
            &= \sum_kh_k2kH_{k-1}(x) \\
            &= \sum_k2(k+1)h_{k+1}H_k(x)
        \end{aligned}
    $$
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    
    See also
    --------
    - in standard monomial basis: [`polyder`][poly.functional.polyder]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermder`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermder.html)
    """
    for _ in range(n):
        h = vechadamard(count(2, 2), islice(h, 1, None))
    return h

def hermantider(h, n=1):
    r"""Return the `n`-th antiderivative of Hermite polynomial series `h`.
    
    $$
        h^{(-n)}
    $$
    
    The first coefficient $h_0$ is always zero. To add an integration constant
    like [`numpy.polynomial.hermite.hermint`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermint.html)
    that would correspond to a lower integration bound
    of $0$ ($\int_0^xh(t)\,\mathrm{d}t$) resulting in $h^{(-n)}(0)=0$
    subtract $h^{(-n)}(0)$.
    
    Notes
    -----
    $$
        \begin{aligned}
            H_k' &= 2kH_{k-1} \\
            \Rightarrow \quad H_k^{(-1)} &= \frac{H_{k+1}}{2(k+1)} \\
            \int h(x)\,\mathrm{d}x &= \int\sum_kh_kH_k(x)\,\mathrm{d}x \\
            &= \sum_kh_k\int H_k(x)\,\mathrm{d}x \\
            &= \sum_kh_k\frac{H_{k+1}(x)}{2(k+1)} \\
            &= \sum_k\frac{h_{k-1}}{2k}H_k(x)
        \end{aligned}
    $$
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    
    See also
    --------
    - in standard monomial basis: [`polyantider`][poly.functional.polyantider]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermint`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermint.html)
    """
    for _ in range(n):
        h = vecrshift(vechadamardtruediv(h, count(2, 2)), 1)
    return h


#sympy
def hermsympify(h, symbol=spx):
    """Return the Hermite polynomial series `h` as a [`sympy.Poly`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly)."""
    return sum((hi*sp.hermite_poly(i, symbol, polys=True) for i, hi in enumerate(h)), start=sp.Poly(0, symbol))
