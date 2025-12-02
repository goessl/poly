from math import comb
from operator import neg, mul
from itertools import repeat, chain, islice, count, starmap, tee, zip_longest
from functools import reduce, cache
from fractions import Fraction
from ..standard import polyzero, polyval, polyadd, polyaddc, polysub, polyscalarmul, polymulx
from .hilbert_space import hermweights, hermweightis
from .vector_space import hermadd, hermaddc, hermscalarmul
from vector import veczero, vecbasis, vecrand, vecrandn



__all__ = (#creation
           'H0', 'H1', 'H2',
           'herm', 'herm_recursive', 'herm_iterative', 'herm_explicit', 'herms',
           'hermzero', 'hermone', 'hermx',
           'hermmono', 'hermmono_recursive', 'hermmono_iterative', 'hermmono_explicit', 'hermmonos',
           'hermrand', 'hermrandn', 'hermfromroots',
           #conversion
           'herm2poly', 'herm2poly_naive', 'herm2poly_iterative', 'herm2poly_clenshaw',
           'poly2herm', 'poly2herm_naive', 'poly2herm_iterative', 'poly2herm_horner',
           #evaluation
           'hermval', 'hermval_naive', 'hermval_iterative', 'hermval_clenshaw',
           'hermvals', 'hermvalzeros', 'hermvalzero',
           #arithmetic
           'hermmul', 'hermmul_naive', 'hermmulx', 'hermmuln',
           'hermpow', 'hermpow_naive', 'hermpow_binary', 'hermmulpow')



#creation
H0 = (1,)
r"""Zero-th Hermite polynomial in standard monomial basis.

$$
    H_0(x)=1 \qquad H_0^{(P)}=\begin{pmatrix} 1 \end{pmatrix}
$$

A tuple with a single one: `(1,)`.
"""

H1 = (0, 2)
r"""First Hermite polynomial in standard monomial basis.

$$
    H_1(x)=2x \qquad H_1^{(P)}=\begin{pmatrix} 0 \\ 2 \end{pmatrix}
$$

A tuple with a zero and a two: `(0, 2)`.
"""

H2 = (-2, 0, 4)
r"""Second Hermite polynomial in standard monomial basis.

$$
    H_2(x)=4x^2-2 \qquad H_2^{(P)}=\begin{pmatrix} -2 \\ 0 \\ 4 \end{pmatrix}
$$

A tuple with a negative two, a zero and a four: `(-2, 0, 4)`.
"""

def herm(n, method='iterative'):
    r"""Return the `n`-th Hermite polynomial in standard monomial basis.
    
    $$
        H_n^{(P)}
    $$
    
    Available methods are
    
    - [`recursive`][poly.hermite.herm_recursive],
    - [`iterative`][poly.hermite.herm_iterative] &
    - [`explicit`][poly.hermite.herm_explicit].
    
    Notes
    -----
    [`numpy.polynomial.hermite.herm2poly(vecbasis(n))`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.herm2poly.html) differs at $n\geq29$.
    
    See also
    --------
    - implementations: [`herm_recursive`][poly.hermite.herm_recursive],
    [`herm_iterative`][poly.hermite.herm_iterative],
    [`herm_explicit`][poly.hermite.herm_explicit]
    - for all Hermite polynomials: [`herms`][poly.hermite.herms]
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
    - for any implementation: [`herm`][poly.hermite.herm]
    - other implementations: [`herm_iterative`][poly.hermite.herm_iterative],
    [`herm_explicit`][poly.hermite.herm_explicit]
    
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
    - for any implementation: [`herm`][poly.hermite.herm]
    - other implementations: [`herm_recursive`][poly.hermite.herm_recursive],
    [`herm_explicit`][poly.hermite.herm_explicit]
    - uses: [`herms`][poly.hermite.herms]
    
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
    - for any implementation: [`herm`][poly.hermite.herm]
    - other implementations: [`herm_recursive`][poly.hermite.herm_recursive],
    [`herm_iterative`][poly.hermite.herm_iterative]
    
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
    (H_0^{(P)})_{n\in\mathbb{N}_0} = (H_0^{(P)}, H_2^{(P)}, H_2^{(P)}, \dots)
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
- in standard monomial basis: [`polyzero`][poly.standard.polyzero]
- wraps: [`vector.veczero`](https://goessl.github.io/vector/functional/#vector.functional.veczero)

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
- in standard monomial basis: [`polyone`][poly.standard.polyone]

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
- in standard monomial basis: [`polyx`][poly.standard.polyx]

References
----------
- `numpy` equivalent: [`numpy.polynomial.hermite.hermx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermx.html)
"""

def hermmono(n, method='explicit'):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        x^n
    $$
    
    The result is a tuple of `Fraction`s.
    
    Available methods are:
    
    - [`recursive`][poly.hermite.hermmono_recursive],
    - [`iterative`][poly.hermite.hermmono_iterative] &
    - [`explicit`][poly.hermite.hermmono_explicit].
    
    See also
    --------
    - implementations: [`hermmono_recursive`][poly.hermite.hermmono_recursive],
    [`hermmono_iterative`][poly.hermite.hermmono_iterative],
    [`hermmono_explicit`][poly.hermite.hermmono_explicit]
    - for all monomials: [`hermmonos`][poly.hermite.hermmonos]
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
    
    Uses [`hermmulx`][poly.hermite.hermmulx] recursively.
    
    Cached, therefore fast if used repeatedly.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite.hermmono]
    - other implementations: [`hermmono_iterative`][poly.hermite.hermmono_iterative],
    [`hermmono_explicit`][poly.hermite.hermmono_explicit]
    """
    if n == 0:
        return (Fraction(1),) #hermone but as `Fraction`
    return hermmulx(hermmono_recursive(n-1))

@cache
def hermmono_iterative(n):
    """Return `x^n` as Hermite polynomial series.
    
    $$
        x^n
    $$
    
    Uses [`hermmulx`][poly.hermite.hermmulx] iteratively.
    
    The result is a tuple of `Fraction`s.
    
    See also
    --------
    - for any implementation: [`hermmono`][poly.hermite.hermmono]
    - other implementations: [`hermmono_recursive`][poly.hermite.hermmono_recursive],
    [`hermmono_explicit`][poly.hermite.hermmono_explicit]
    - uses: [`hermmonos`][poly.hermite.hermmonos]
    """
    return next(islice(hermmonos(), n, None))

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
    - for any implementation: [`hermmono`][poly.hermite.hermmono]
    - other implementations: [`hermmono_recursive`][poly.hermite.hermmono_recursive],
    [`hermmono_explicit`][poly.hermite.hermmono_explicit]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Inverse explicit expression](https://en.wikipedia.org/wiki/Hermite_polynomials#Inverse_explicit_expression)
    """
    r = [Fraction()] * (n + 1)
    r[n] = Fraction(1, 2**n)
    for m in range(1, n//2+1):
        r[n-2*m] = r[n-2*(m-1)] * ((n-2*m+1)*(n-2*m+2)) / m
    return tuple(r)

def hermmonos():
    r"""Yield standard monomials `x^n` as Hermite polynomial series.
    
    $$
        \left(x^n\right)_{n\in\mathbb{N}_0} = (1, x, x^2, \dots)
    $$
    
    Uses [`hermmulx`][poly.hermite.hermmulx] iteratively.
    
    See also
    --------
    - used in: [`hermmono_iterative`][poly.hermite.hermmono_iterative]
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
    - in standard monomial basis: [`polyrand`][poly.standard.polyrand]
    - wraps: [`vector.vecrand`](https://goessl.github.io/vector/functional/#vector.functional.vecrand)
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
    - in standard monomial basis: [`polyrandn`][poly.standard.polyrandn]
    - wraps: [`vector.vecrandn`](https://goessl.github.io/vector/functional/#vector.functional.vecrandn)
    """
    return vecrandn(n+1, normed=normed, mu=mu, sigma=sigma, weights=hermweights())

def hermfromroots(*xs):
    r"""Return the Hermite polynomial series with the given roots.
    
    $$
        \prod_k(x-x_k)
    $$
    
    See also
    --------
    - in standard monomial basis: [`polyfromroots`][poly.standard.polyfromroots]
    """
    return hermmul(*zip(map(neg, xs), repeat(Fraction(1, 2))))


#conversion
def herm2poly(h, method='iterative'):
    """Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Available methods are 
    
    - [`naive`][poly.hermite.herm2poly_naive],
    - [`iterative`][poly.hermite.herm2poly_iterative] &
    - [`clenshaw`][poly.hermite.herm2poly_clenshaw].
    
    See also
    --------
    - implementations: [`herm2poly_naive`][poly.hermite.herm2poly_naive],
    [`herm2poly_iterative`][poly.hermite.herm2poly_iterative],
    [`herm2poly_clenshaw`][poly.hermite.herm2poly_clenshaw]
    
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
    r"""Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Uses naive Hermite polynomial creation by [`herm_explicit`][poly.hermite.herm_explicit].
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $\frac{n(n+1)}{2}$ scalar additions (`add`) &
    - $\frac{(n+1)(n+2)}{2}$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`herm2poly`][poly.hermite.herm2poly]
    - other implementations: [`herm2poly_iterative`][poly.hermite.herm2poly_iterative],
    [`herm2poly_clenshaw`][poly.hermite.herm2poly_clenshaw]
    """
    return polyadd(*(polyscalarmul(hi, herm_explicit(i)) for i, hi in enumerate(h)))

def herm2poly_iterative(h):
    r"""Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Uses iterative Hermite polynomial creation like
    [`herm_iterative`][poly.hermite.herm_iterative].
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $\frac{n(n+1)}{2}$ scalar additions (`add`) &
    - $\frac{(n+1)(n+2)}{2}$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`herm2poly`][poly.hermite.herm2poly]
    - other implementations: [`herm2poly_naive`][poly.hermite.herm2poly_naive],
    [`herm2poly_clenshaw`][poly.hermite.herm2poly_clenshaw]
    - uses: [`herms`][poly.hermite.herms]
    """
    return polyadd(*starmap(polyscalarmul, zip(h, herms())))

def herm2poly_clenshaw(h):
    r"""Return a Hermite polynomial series in standard monomial basis.
    
    $$
        h^{(P)}
    $$
    
    Uses the Clenshaw algorithm like
    [`hermval_clenshaw`][poly.hermite.hermval_clenshaw].
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar additions (`add`),
    - $\begin{cases}\frac{n(n-1)}{2} & n\geq0 \\ 0 & n\leq0\end{cases}$ scalar subtractions (`sub`) &
    - $\begin{cases}n^2 & n\geq0 \\ 0 & n\leq0\end{cases}$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`herm2poly`][poly.hermite.herm2poly]
    - other implementations [`herm2poly_naive`][poly.hermite.herm2poly_naive],
    [`herm2poly_iterative`][poly.hermite.herm2poly_iterative]
    """
    h = tuple(h)
    if not h:
        return polyzero
    a, b = polyzero, polyzero
    for n, hn in reversed(tuple(enumerate(h[1:], 1))):
        a, b = polyaddc(polysub(polyscalarmul(2, polymulx(a)), polyscalarmul(2*(n+1), b)), hn), a
    return polyaddc(polysub(polyscalarmul(2, polymulx(a)), polyscalarmul(2, b)), h[0])

def poly2herm(p, method='iterative'):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Available methods are 
    
    - [`naive`][poly.hermite.poly2herm_naive],
    - [`iterative`][poly.hermite.poly2herm_iterative] &
    - [`horner`][poly.hermite.poly2herm_horner].
    
    See also
    --------
    - implementations: [`poly2herm_naive`][poly.hermite.poly2herm_naive],
    [`poly2herm_iterative`][poly.hermite.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite.poly2herm_horner]
    
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
    method of [`hermmono`][poly.hermite.hermmono].
    
    See also
    --------
    - any implementation: [`poly2herm`][poly.hermite.poly2herm]
    - other implementations: [`poly2herm_iterative`][poly.hermite.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite.poly2herm_horner]
    """
    return hermadd(*(hermscalarmul(pi, hermmono(i)) for i, pi in enumerate(p)))

def poly2herm_iterative(p):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Uses iterative monomial creation of [`hermmonos`][poly.hermite.hermmonos].
    
    See also
    --------
    - for any implementation: [`poly2herm`][poly.hermite.poly2herm]
    - other implementations: [`poly2herm_iterative`][poly.hermite.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite.poly2herm_horner]
    """
    return hermadd(*starmap(hermscalarmul, zip(p, hermmonos())))

def poly2herm_horner(p):
    """Return a standard monomial basis polynomials as a Hermite polynomial series.
    
    $$
        p^{(H)}
    $$
    
    Uses Horner's method.
    
    `p` must be reversible.
    
    See also
    --------
    - for any implementation: [`poly2herm`][poly.hermite.poly2herm]
    - other implementations: [`poly2herm_iterative`][poly.hermite.poly2herm_iterative],
    [`poly2herm_horner`][poly.hermite.poly2herm_horner]
    
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
    
    - [`naive`][poly.hermite.hermval_naive],
    - [`iterative`][poly.hermite.hermval_iterative] &
    - [`clenshaw`][poly.hermite.hermval_clenshaw] (`h` must be reversible).
    
    See also
    --------
    - implementations: [`hermval_naive`][poly.hermite.hermval_naive],
    [`hermval_iterative`][poly.hermite.hermval_iterative],
    [`hermval_clenshaw`][poly.hermite.hermval_clenshaw]
    - in standard monomial basis: [`polyval`][poly.standard.polyval]
    
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
    
    Converts $h$ to standard monomial basis and evaluates with [`polyval`][poly.standard.polyval].
    
    See also
    --------
    - for any implementation: [`hermval`][poly.hermite.hermval]
    - other implementations: [`hermval_iterative`][poly.hermite.hermval_iterative],
    [`hermval_clenshaw`][poly.hermite.hermval_clenshaw]
    - in standard monomial basis: [`polyval_naive`][poly.standard.polyval_naive]
    """
    return polyval(herm2poly(h), x)

def hermval_iterative(h, x):
    """Return the value of Hermite polynomial series `h` evaluated at point `x`.
    
    $$
        h(x)
    $$
    
    Uses iterative $H_n(x)$ generation.
    
    See also
    --------
    - for any implementation: [`hermval`][poly.hermite.hermval]
    - other implementations: [`hermval_naive`][poly.hermite.hermval_naive],
    [`hermval_clenshaw`][poly.hermite.hermval_clenshaw]
    - uses: [`hermvals`][poly.hermite.hermvals]
    - in standard monomial basis: [`polyval_iterative`][poly.standard.polyval_iterative]
    """
    return sum(map(mul, h, hermvals(x)), start=type(x)(0))

def hermval_clenshaw(h, x):
    """Return the value of Hermite polynomial series `h` evaluated at point `x`.
    
    $$
        h(x)
    $$
    
    Uses the Clenshaw algorithm.
    
    `h` must be reversible.
    
    See also
    --------
    - for any implementation: [`hermval`][poly.hermite.hermval]
    - other implementations: [`hermval_naive`][poly.hermite.hermval_naive],
    [`hermval_iterative`][poly.hermite.hermval_iterative]
    - in standard monomial basis: [`polyval_horner`][poly.standard.polyval_horner]
    
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

def hermvals(x):
    r"""Yield the values of Hermite polynomials evaluated at point `x`.
    
    $$
        (H_n(x))_{n\in\mathbb{N}_0} = (H_0(x), H_1(x), H_2(x), \dots)
    $$
    
    Uses the recurrence relation $H_n(x)=2xH_{n-1}(x)-2(n-1)H_{n-2}(x)$ iteratively.
    
    See also
    --------
    - used in: [`hermval_iterative`][poly.hermite.hermval_iterative]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    yield (Hkm1 := type(x)(1))
    yield (Hk := 2*x)
    for k in count(1):
        Hkm1, Hk = Hk, 2*(Hk*x - k*Hkm1)
        yield Hk

def hermvalzeros():
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
    - for any series: [`hermvalzero`][poly.hermite.hermvalzero]
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
    - for any point: [`hermval`][poly.hermite.hermval]
    """
    return sum(starmap(mul, islice(zip(h, hermvalzeros()), 0, None, 2)))


#arithmetic
def hermmul(*hs, method='naive'):
    r"""Return the product of Hermite polynomial series.
    
    $$
        h_0 h_1 \cdots
    $$
    
    Available methods are
    
    - [`naive`][poly.hermite.hermmul_naive].
    
    TODO: Karatsuba?
    
    See also
    --------
    - implementations: [`hermmul_naive`][poly.hermite.hermmul_naive]
    - for monomial factor: [`hermmulx`][poly.hermite.hermmulx]
    - for Hermite polynomial factor: [`hermmulx`][poly.hermite.hermmuln]
    - in standard monomial basis: [`polymul`][poly.standard.polymul]
    
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
    
    TODO: make g only iterable (dynamicly extending list)
    
    Complexity
    ----------
    For two Hermite polynomial series of degrees $n$ & $m$ (let $N=\min\{n, m\}$, $M=\max\{n, m\}$) there will be
    
    - $\begin{cases}\frac{(N+1)(N+2)(3M+3-N)}{6}-(n+m+1) & n\geq0\land m\geq0 \\ 0 & n<0\lor m<0 \end{cases}$ scalar additions (`add`) &
    - $\begin{cases}\frac{(N+1)(N+2)(3M+3-N)}{3} & n\geq0\land m\geq0 \\ 0 & n<0\lor m<0 \end{cases}$ scalar multiplications (`mul`).
    
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
    - for any implementation: [`hermmul`][poly.hermite.hermmul]
    - in standard monomial basis: [`polymul_naive`][poly.standard.polymul_naive]
    
    References
    ----------
    - [1] Zhi-yuan Huang & Jia-an Yan: Introduction to Infinite Dimensional Stochastic Analysis. [10.1007/978-94-011-4108-6](https://doi.org/10.1007/978-94-011-4108-6)
    - [Wikipedia - Hermite polynomials - Definition](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition)
    """
    if not g or not h:
        return hermzero
    sentinel = object()
    r = [sentinel] * (len(g)+len(h)-1)
    for i, gi in enumerate(g):
        for j, hj in enumerate(h):
            for k, f in enumerate(islice(hermweightis(), min(i, j)+1)):
                if r[i+j-2*k] is sentinel:
                    r[i+j-2*k] = f * comb(i, k)*comb(j, k) * gi*hj
                else:
                    r[i+j-2*k] += f * comb(i, k)*comb(j, k) * gi*hj
    return tuple(r)

def hermmulx(h):
    r"""Return the product of Hermite polynomial series `h` and a monomial of degree `n`.
    
    $$
        h(x)x^n
    $$
    
    More efficient than `hermmul(h, hermx)`.
    
    TODO: Type safety & complexity
    
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
    - for Hermite polynomial series as factor: [`hermmul`][poly.hermite.hermmul]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermmulx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermmulx.html)
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    im1, ip1 = tee(h)
    _, h1 = next(ip1, None), next(ip1, 0)
    return tuple(chain((h1,), (hkm1/2+(k+1)*hkp1 for k, (hkm1, hkp1) in
            enumerate(zip_longest(im1, ip1, fillvalue=0), 1))))

def hermmuln(h, n):
    r"""Return the product of Hermite polynomial series `h` and the `n`-th Hermite polynomial.
    
    $$
        hH_n
    $$
    
    More efficient than `hermmul(h, vecbasis(n))`.
    
    TODO: Type safety & complexity
    
    See also
    --------
    - for a Hermite polynomial series factor: [`hermmul`][poly.hermite.hermmul]
    """
    r = []
    for i, hi in enumerate(h):
        for k, f in enumerate(islice(hermweightis(), min(i, n)+1)):
            r.extend([0] * (i+n-2*k - len(r) + 1))
            r[i+n-2*k] += f * comb(i, k)*comb(n, k) * hi
    return tuple(r)

def hermpow(h, n, method='binary'):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    Available methods are 
    
    - [`naive`][poly.hermite.hermpow_naive] &
    - [`binary`][poly.hermite.hermpow_binary].
    
    See also
    --------
    - implementations: [`hermpow_naive`][poly.hermite.hermpow_naive],
    [`hermpow_binary`][poly.hermite.hermpow_binary],
    - in standard monomial basis: [`polypow`][poly.standard.polypow]
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
    - for any implementation: [`hermpow`][poly.hermite.hermpow]
    - other implementations: [`hermpow_binary`][poly.hermite.hermpow_binary]
    - in standard monomial basis: [`polypow_naive`][poly.standard.polypow_naive]
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
    - for any implementations: [`hermpow`][poly.hermite.hermpow]
    - other implementations[`hermpow_naive`][poly.hermite.hermpow_naive]
    - in standard monomial basis: [`polypow_binary`][poly.standard.polypow_binary]
    
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
