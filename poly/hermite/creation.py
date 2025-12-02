from operator import neg
from itertools import count, islice, repeat
from functools import cache
from fractions import Fraction
from ..standard import polysub, polyscalarmul, polymulx
from .hilbert_space import hermweights
from .arithmetic import hermmul, hermmulx
from vector import veczero, vecrand, vecrandn



__all__ = ('H0', 'H1', 'H2',
           'herm', 'herm_recursive', 'herm_iterative', 'herm_explicit', 'herms',
           'hermzero', 'hermone', 'hermx',
           'hermmono', 'hermmono_recursive', 'hermmono_iterative', 'hermmono_explicit', 'hermmonos',
           'hermrand', 'hermrandn', 'hermfromroots')



H0 = (1,)
r"""Zero-th Hermite polynomial in standard monomial basis.

$$
    H_0(x)=1 \qquad H_0^{(P)}=\begin{pmatrix} 1 \end{pmatrix}
$$

A tuple with a single one: `(1,)`.

See also
--------
- other constants: [`H1`][poly.hermite.creation.H1], [`H2`][poly.hermite.creation.H2]
- for any degree: [`herm`][poly.hermite.creation.herm]
"""

H1 = (0, 2)
r"""First Hermite polynomial in standard monomial basis.

$$
    H_1(x)=2x \qquad H_1^{(P)}=\begin{pmatrix} 0 \\ 2 \end{pmatrix}
$$

A tuple with a zero and a two: `(0, 2)`.

See also
--------
- other constants: [`H0`][poly.hermite.creation.H0], [`H2`][poly.hermite.creation.H2]
- for any degree: [`herm`][poly.hermite.creation.herm]
"""

H2 = (-2, 0, 4)
r"""Second Hermite polynomial in standard monomial basis.

$$
    H_2(x)=4x^2-2 \qquad H_2^{(P)}=\begin{pmatrix} -2 \\ 0 \\ 4 \end{pmatrix}
$$

A tuple with a negative two, a zero and a four: `(-2, 0, 4)`.

See also
--------
- other constants: [`H0`][poly.hermite.creation.H0], [`H1`][poly.hermite.creation.H1]
- for any degree: [`herm`][poly.hermite.creation.herm]
"""

def herm(n, method='iterative'):
    r"""Return the `n`-th Hermite polynomial in standard monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Available methods are
    
    - [`recursive`][poly.hermite.creation.herm_recursive],
    - [`iterative`][poly.hermite.creation.herm_iterative] &
    - [`explicit`][poly.hermite.creation.herm_explicit].
    
    Notes
    -----
    [`numpy.polynomial.hermite.herm2poly(vecbasis(n))`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.herm2poly.html) differs at $n\geq29$.
    
    See also
    --------
    - constants: [`H0`][poly.hermite.creation.H0], [`H1`][poly.hermite.creation.H1], [`H2`][poly.hermite.creation.H2]
    - implementations: [`herm_recursive`][poly.hermite.creation.herm_recursive],
    [`herm_iterative`][poly.hermite.creation.herm_iterative],
    [`herm_explicit`][poly.hermite.creation.herm_explicit]
    - for all Hermite polynomials: [`herms`][poly.hermite.creation.herms]
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
    """Return the `n`-th Hermite polynomial in standard  monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ recursively.
    
    Cached, therefore fast if used repeatedly.
    
    See also
    --------
    - for any implementation: [`herm`][poly.hermite.creation.herm]
    - other implementations: [`herm_iterative`][poly.hermite.creation.herm_iterative],
    [`herm_explicit`][poly.hermite.creation.herm_explicit]
    
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
    """Return the `n`-th Hermite polynomial in standard  monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ iteratively.
    
    See also
    --------
    - for any implementation: [`herm`][poly.hermite.creation.herm]
    - other implementations: [`herm_recursive`][poly.hermite.creation.herm_recursive],
    [`herm_explicit`][poly.hermite.creation.herm_explicit]
    - uses: [`herms`][poly.hermite.creation.herms]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    return next(islice(herms(), n, None))

def herm_explicit(n):
    r"""Return the `n`-th Hermite polynomial in standard  monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Uses the explicit expression $H_n(x)=n!\sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}\frac{(-1)^m}{m!(n-2m)!}(2x)^{n-2m}$.
    
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
    - for any implementation: [`herm`][poly.hermite.creation.herm]
    - other implementations: [`herm_recursive`][poly.hermite.creation.herm_recursive],
    [`herm_iterative`][poly.hermite.creation.herm_iterative]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Explicit expression](https://en.wikipedia.org/wiki/Hermite_polynomials#Explicit_expression)
    """
    r = [0] * (n+1)
    r[-1] = 2**n
    for m in range(1, n//2+1):
        r[n-2*m] = -r[n-2*(m-1)] * (n-2*m+1)*(n-2*m+2) // (4*m)
    return tuple(r)

def herms():
    r"""Yield the Hermite polynomials in standard monomial basis.
    
    $$
        (H_n^{(P)})_{n\in\mathbb{N}_0} = (H_0^{(P)}, H_2^{(P)}, H_2^{(P)}, \dots)
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ iteratively.
    
    See also
    --------
    - used in: [`herm_iterative`][poly.hermite.herm_iterative]
    
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
r"""Zero Hermite polynomial series.

$$
    0 \qquad 0^{(H)}=\vec{0}
$$

An empty tuple: `()`.

See also
--------
- other constants: [`hermone`][poly.hermite.creation.hermone], [`hermx`][poly.hermite.creation.hermx]
- for any monomial: [`hermmono`][poly.hermite.creation.hermmono]
- wraps: [`vector.veczero`](https://goessl.github.io/vector/functional/#vector.functional.creation.veczero)

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermzero`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermzero.html)
"""

hermone = (1, )
r"""Constant one Hermite polynomial series.

$$
    1=H_0 \qquad 1^{(H)}=\begin{pmatrix} 1 \end{pmatrix}
$$

A tuple with a single one: `(1,)`.

Notes
-----
A `Fraction` `tuple` (`(Fraction(1), )`) would be more consitent with [`hermx`][poly.hermite.hermx], but would then require `Fraction` arithmetic for multiplicative functions, ([`hermmul`][poly.hermite.hermmul], [`hermpow`][poly.hermite.hermpow]).

See also
--------
- other constants: [`hermzero`][poly.hermite.creation.hermzero], [`hermx`][poly.hermite.creation.hermx]
- for any monomial: [`hermmono`][poly.hermite.creation.hermmono]

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermone`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermone.html)
"""

hermx = (Fraction(0), Fraction(1, 2))
r"""Identity Hermite polynomial series.

$$
    x=\frac{1}{2}H_1(x) \qquad x^{(H)}=\begin{pmatrix} 0 \\ \frac{1}{2} \end{pmatrix}
$$

A tuple with a zero and a half: `(0, 1/2)`.

See also
--------
- other constants: [`hermzero`][poly.hermite.creation.hermzero], [`hermone`][poly.hermite.creation.hermone]
- for any monomial: [`hermmono`][poly.hermite.creation.hermmono]

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermx.html)
"""

def hermmono(n, method='explicit'):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        (x^n)^{(H)}
    $$
    
    The result is a tuple of `Fraction`s.
    
    Available methods are:
    
    - [`recursive`][poly.hermite.creation.hermmono_recursive],
    - [`iterative`][poly.hermite.creation.hermmono_iterative] &
    - [`explicit`][poly.hermite.creation.hermmono_explicit].
    
    See also
    --------
    - constants: [`hermzero`][poly.hermite.creation.hermzero], [`hermone`][poly.hermite.creation.hermone], [`hermx`][poly.hermite.creation.hermx]
    - implementations: [`hermmono_recursive`][poly.hermite.creation.hermmono_recursive],
    [`hermmono_iterative`][poly.hermite.creation.hermmono_iterative],
    [`hermmono_explicit`][poly.hermite.creation.hermmono_explicit]
    - for all monomials: [`hermmonos`][poly.hermite.creation.hermmonos]
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
        x^n
    $$
    
    Uses [`hermmulx`][poly.hermite.arithmetic.hermmulx] recursively.
    
    Cached, therefore fast if used repeatedly.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite.creation.hermmono]
    - other implementations: [`hermmono_iterative`][poly.hermite.creation.hermmono_iterative],
    [`hermmono_explicit`][poly.hermite.creation.hermmono_explicit]
    """
    if n == 0:
        return (Fraction(1),) #hermone but as `Fraction`
    return hermmulx(hermmono_recursive(n-1))

def hermmono_iterative(n):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        x^n
    $$
    
    Uses [`hermmulx`][poly.hermite.arithmetic.hermmulx] iteratively.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite.creation.hermmono]
    - other implementations: [`hermmono_recursive`][poly.hermite.creation.hermmono_recursive],
    [`hermmono_explicit`][poly.hermite.creation.hermmono_explicit]
    - uses: [`hermmonos`][poly.hermite.creation.hermmonos]
    """
    return next(hermmonos(n))

def hermmono_explicit(n):
    r"""Return `x^n` as Hermite polynomial series.
    
    $$
        x^n
    $$
    
    Uses the explicit expression $x^n=\frac{n!}{2^n}\sum_{m=0}^{\lfloor\frac{n}{2}\rfloor}\frac{1}{m!(n-2m)!}H_{n-2m}(x)$.
    
    The result is a tuple of `Fraction`s.
    
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
    - for any implementation: [`hermmono`][poly.hermite.creation.hermmono]
    - other implementations: [`hermmono_recursive`][poly.hermite.creation.hermmono_recursive],
    [`hermmono_explicit`][poly.hermite.creation.hermmono_explicit]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Inverse explicit expression](https://en.wikipedia.org/wiki/Hermite_polynomials#Inverse_explicit_expression)
    """
    r = [Fraction()] * (n + 1)
    r[n] = Fraction(1, 2**n)
    for m in range(1, n//2+1):
        r[n-2*m] = r[n-2*(m-1)] * ((n-2*m+1)*(n-2*m+2)) / m
    return tuple(r)

def hermmonos(start=0):
    r"""Yield standard monomials `x^n` as Hermite polynomial series.
    
    $$
        \left(x^n\right)_{n\in\mathbb{N}_0} = (1, x, x^2, \dots)
    $$
    
    Uses [`hermmulx`][poly.hermite.arithmetic.hermmulx] iteratively.
    
    See also
    --------
    - used in: [`hermmono_iterative`][poly.hermite.creation.hermmono_iterative]
    """
    h = (Fraction(1),) #hermone but as `Fraction`
    for _ in range(start):
        h = hermmulx(h)
    yield h
    for _ in count(start):
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
    - wraps: [`vector.vecrand`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecrand)
    """
    return vecrand(n+1)

def hermrandn(n, normed=True, mu=0, sigma=1):
    r"""Return a random Hermite polynomial series of degree `n`.
    
    $$
        \sum_{k=0}^na_kH_k(x) \quad a_k\sim\mathcal{N}(\mu, \sigma)
    $$
    
    The coefficients are sampled from a normal distribution.
    
    Normed with respect to the Hermite polynomial norm $\int_\mathbb{R}h^2(x)e^{-x^2}\,\mathrm{d}x$.
    
    See also
    --------
    - wraps: [`vector.vecrandn`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecrandn)
    """
    return vecrandn(n+1, normed=normed, mu=mu, sigma=sigma, weights=hermweights())

def hermfromroots(*xs):
    r"""Return the Hermite polynomial series with the given roots.
    
    $$
        \prod_k(x-x_k)
    $$
    """
    return hermmul(*zip(map(neg, xs), repeat(Fraction(1, 2))))
