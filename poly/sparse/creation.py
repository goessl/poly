from vector.sparse.creation import vecszero, vecsbasis, vecsbases, vecsrand, vecsrandn
from .arithmetic import polyssub, polysscalarmul, polysmulx



__all__ = ('polyszero', 'polysone', 'polysx', 'polysmono', 'polysmonos', 'polysrand', 'polysrandn', 'polysfromroots')



polyszero = vecszero
"""Zero polynomial.

$$
    0
$$

An empty `dict`: `{}`.

See also
--------
- other constants: [`polysone`][poly.sparse.creation.polysone], [`polysx`][poly.sparse.creation.polysx]
- for any degree: [`polysmono`][poly.sparse.creation.polysmono]
- wraps: [`vector.vecszero`](https://goessl.github.io/vector/sparse/#vector.sparse.creation.vecszero)
"""

polysone = {0:1}
"""Constant one polynomial.

$$
    1
$$

See also
--------
- other constants: [`polyszero`][poly.sparse.creation.polyszero], [`polysx`][poly.sparse.creation.polysx]
- for any degree: [`polysmono`][poly.sparse.creation.polysmono]
"""

polysx = {1:1}
"""Identity polynomial.

$$
    x
$$

See also
--------
- other constants: [`polyszero`][poly.sparse.creation.polyszero], [`polysone`][poly.sparse.creation.polysone]
- for any degree: [`polysmono`][poly.sparse.creation.polysmono]
"""

def polysmono(i, c=1):
    r"""Return a monomial of degree `n`.
    
    $$
        cx^n
    $$
    
    See also
    --------
    - other constants: [`polyszero`][poly.sparse.creation.polyszero], [`polysone`][poly.sparse.creation.polysone], [`polysx`][poly.sparse.creation.polysx]
    - for all monomials: [`polysmonos`][poly.sparse.creation.polysmonos]
    - wraps: [`vector.vecsbasis`](https://goessl.github.io/vector/sparse/#vector.sparse.creation.vecsbasis)
    """
    return vecsbasis(i, c=c) if i>=0 else polyszero

def polysmonos(start=0, c=1):
    r"""Yield all monomials.
    
    $$
        \left(x^n\right)_\mathbb{n\in\mathbb{N_0}} = \left(1, x, x^2, \dots \right)
    $$
    
    See also
    --------
    - for single monomial: [`polysmono`][poly.sparse.creation.polysmono]
    - wraps: [`vector.vecsbases`](https://goessl.github.io/vector/sparse/#vector.sparse.creation.vecsbases)
    """
    yield from vecsbases(start=start, c=c)

def polysrand(n):
    r"""Return a random polynomial of degree `n`.
    
    $$
        \sum_ka_kx^k \quad a_k\sim\mathcal{U}^n([0, 1[)
    $$
    
    The coefficients are sampled from a uniform distribution in `[0, 1[`.
    
    See also
    --------
    - wraps: [`vector.vecsrand`](https://goessl.github.io/vector/sparse/#vector.sparse.creation.vecsrand)
    """
    return vecsrand(n+1) if n>=0 else polyszero

def polysrandn(n, normed=True, mu=0, sigma=1, weights=None):
    r"""Return a random polynomial of degree `n`.
    
    $$
        \sum_{k=0}^na_kx^k \quad a_k\sim\mathcal{N}(\mu, \sigma)
    $$
    
    The coefficients are sampled from a normal distribution.
    
    Normed with respect to the euclidian vector norm $\sum_k|a_k|^2$.
    
    See also
    --------
    - wraps: [`vector.vecsrandn`](https://goessl.github.io/vector/functional/#vector.functional.creation.vecsrandn)
    """
    return vecsrandn(n+1, normed=normed, mu=mu, sigma=sigma) if n>=0 else polyszero

def polysfromroots(*xs, one=1):
    r"""Return the polynomial with the given roots.
    
    $$
        \prod_k(x-x_k)
    $$
    """
    r = {0:one}
    for x in xs:
        r = polyssub(polysmulx(r), polysscalarmul(x, r))
    return r
