from itertools import count, islice
from ..standard import polyval
from .conversion import herm2poly
from vector import vecdot



__all__ = ('hermval', 'hermval_naive', 'hermval_iterative', 'hermval_clenshaw',
           'hermvals', 'hermvalzeros', 'hermvalzero')



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
    return vecdot(h, hermvals(x), zero=type(x)(0))

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

def hermvalzero(h, zero=0):
    """Return the value of Hermite polynomial series `h` evaluated at point 0.
    
    $$
        h(0)
    $$
    
    More efficient than `hermval(h, 0)`.
    
    See also
    --------
    - for any point: [`hermval`][poly.hermite.hermval]
    """
    return vecdot(islice(h, 0, None, 2), islice(hermvalzeros(), 0, None, 2), zero=zero)
