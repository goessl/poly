from .arithmetic import polysub, polyscalarmul, polymulx
from vector import vecbasis, vecbases, vecrand, vecrandn



__all__ = ('polyzero', 'polyone', 'polyx', 'polymono', 'polymonos',
           'polyrand', 'polyrandn', 'polyfromroots')



polyzero = ()
r"""Zero polynomial.

$$
    0 \qquad 0^{(P)}=\vec{0}
$$

An empty tuple: `()`.

Notes
-----
Why give the zero polynomial a distinguished representation (not just `[0]` like [`numpy.polynomial`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyzero.html))?

The additive neutral element seems worth handling exceptionally.
It is mathematically different (different degree)
and results in functions like [`polymul`][poly.standard.arithmetic.polymul] being more time and memory efficient.

See also
--------
- other constants: [`polyone`][poly.standard.polyone], [`polyx`][poly.standard.polyx]
- for any degree: [`polymono`][poly.standard.polymono]
- wraps: [`vector.veczero`](https://goessl.github.io/vector/functional/#vector.functional.veczero)

References
----------
- `numpy` equivalent: [`numpy.polynomial.polynomial.polyzero`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyzero.html)
"""

polyone = (1, )
r"""Constant one polynomial.

$$
    1 \qquad 1^{(P)}=\begin{pmatrix} 1 \end{pmatrix}
$$

A tuple with a single one: `(1,)`.

See also
--------
- other constants: [`polyzero`][poly.standard.polyzero], [`polyx`][poly.standard.polyx]
- for any degree: [`polymono`][poly.standard.polymono]

References
----------
- `numpy` equivalent: [`numpy.polynomial.polynomial.polyone`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyone.html)
"""

polyx = (0, 1)
r"""Identity polynomial.

$$
    x \qquad x^{(P)}=\begin{pmatrix} 0 \\ 1 \end{pmatrix}
$$

A tuple with a zero and a one `(0, 1)`.

See also
--------
- other constants: [`polyzero`][poly.standard.polyzero], [`polyone`][poly.standard.polyone]
- for any degree: [`polymono`][poly.standard.polymono]

References
----------
- `numpy` equivalent: [`numpy.polynomial.polynomial.polyx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyx.html)
"""

def polymono(n, c=1, zero=0):
    r"""Return a monomial of degree `n`.
    
    $$
        cx^n \qquad (cx^n)^{(P)}=c\vec{e}_n
    $$
    
    Returns a tuple with `n` zeros followed by `c` or [`polyzero`][poly.standard.polyzero] if $n<0$.
    
    See also
    --------
    - constants: [`polyzero`][poly.standard.polyzero], [`polyone`][poly.standard.polyone], [`polyx`][poly.standard.polyx]
    - for all monomials: [`polymonos`][poly.standard.polymonos]
    - wraps: [`vector.vecbasis`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecbasis)
    """
    return vecbasis(n, c=c, zero=zero) if n>=0 else polyzero

def polymonos(start=0, c=1, zero=0):
    r"""Yield all monomials.
    
    $$
        \left(x^n\right)_\mathbb{n\in\mathbb{N_0}} = \left(1, x, x^2, \dots \right)
    $$
    
    See also
    --------
    - for single monomial: [`polymono`][poly.standard.polymono]
    - wraps: [`vector.vecbases`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecbases)
    """
    yield from vecbases(start=start, c=c, zero=zero)

def polyrand(n):
    r"""Return a random polynomial of degree `n`.
    
    $$
        \sum_{k=0}^na_kx^k \quad a_k\sim\mathcal{U}([0, 1[)
    $$
    
    The coefficients are sampled from a uniform distribution in `[0, 1[`.
    
    See also
    --------
    - wraps: [`vector.vecrand`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecrand)
    """
    return vecrand(n+1)

def polyrandn(n, normed=True, mu=0, sigma=1):
    r"""Return a random polynomial of degree `n`.
    
    $$
        \sum_{k=0}^na_kx^k \quad a_k\sim\mathcal{N}(\mu, \sigma)
    $$
    
    The coefficients are sampled from a normal distribution.
    
    Normed with respect to the euclidian vector norm $\sum_k|a_k|^2$.
    
    See also
    --------
    - wraps: [`vector.vecrandn`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecrandn)
    """
    return vecrandn(n+1, normed=normed, mu=mu, sigma=sigma)

def polyfromroots(*xs, one=1):
    r"""Return the polynomial with the given roots.
    
    $$
        \prod_k(x-x_k)
    $$
    
    Complexity
    ----------
    For $n$ roots there will be
    
    - $n$ scalar negations (`neg`),
    - $\frac{n(n-1)}{2}$ scalar additions (`add`) &
    - $\begin{cases}(n+2)(n-1)&n\ge1\\0&n\le1\end{cases}$ scalar multiplications (`mul`).
    
    References
    ----------
    - Recipe: [more_itertools.polynomial_from_roots](https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.polynomial_from_roots)
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polyfromroots`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyfromroots.html)
    """
    #return polymul(*zip(map(neg, xs), repeat(one)), method='naive', one=one)
    r = (one, )
    for x in xs:
        r = polysub(polymulx(r), polyscalarmul(x, r))
    return r
