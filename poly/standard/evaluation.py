from itertools import chain, count, islice, repeat
from .arithmetic import polyadd, polyaddc, polyscalarmul, polymul_naive, polypow_naive, polypows
from operationcounter import MISSING, reduce_default, prod_default, sumprod_default
from vector import veclhadamard



__all__ = ('polyval', 'polyval_naive', 'polyval_iterative', 'polyval_horner', 'polyvals', 'polyvalzero',
           'polycom', 'polycom_naive', 'polycom_iterative', 'polycom_horner',
           'polyshift', 'polyscale')



def polyval(p, x, method='horner'):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Available methods are
    
    - [`naive`][poly.standard.polyval_naive],
    - [`iterative`][poly.standard.polyval_iterative] &
    - [`horner`][poly.standard.polyval_horner] (`p` must be reversible).
    
    See also
    --------
    - implementations: [`polyval_naive`][poly.standard.polyval_naive],
    [`polyval_iterative`][poly.standard.polyval_iterative],
    [`polyval_horner`][poly.standard.polyval_horner]
    - for consecutive monomials: [`polyvals`][poly.standard.polyvals]
    - for $x=0$: [`polyvalzero`][poly.standard.polyvalzero]
    - for polynomial arguments: [`polycom`][poly.standard.polycom]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polyval`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyval.html)
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

def polyval_naive(p, x):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Uses naive repeated multiplication to calculate monomials individually.
    
    See also
    --------
    - for any implementation: [`polyval`][poly.standard.polyval]
    - other implementations:
    [`polyval_iterative`][poly.standard.polyval_iterative],
    [`polyval_horner`][poly.standard.polyval_horner]
    - for polynomial arguments: [`polycom_naive`][poly.standard.polycom_naive]
    
    References
    ----------
    - [more_itertools.polynomial_eval](https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.polynomial_eval)
    - [Wikipedia - Polynomial evaluation](https://en.wikipedia.org/wiki/Polynomial_evaluation)
    - [Wikipedia - Horner's method - Efficiency](https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency)
    """
    p = iter(p)
    a0 = next(p, type(x)(0))
    return sumprod_default(p, chain((x,), map(pow, repeat(x), count(start=2))), initial=a0, default=MISSING)

def polyval_iterative(p, x):
    r"""Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Uses iterative multiplication to calculate monomials consecutively.
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $\begin{cases}n&n\ge0\\0&n\le0\end{cases}$ scalar additions (`add`) &
    - $\begin{cases}2n-1&n>0\\0&n\le0\end{cases}$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`polyval`][poly.standard.polyval]
    - other implementations:
    [`polyval_naive`][poly.standard.polyval_naive],
    [`polyval_horner`][poly.standard.polyval_horner]
    - uses: [`polyvals`][poly.standard.polyvals]
    - for polynomial arguments: [`polycom_iterative`][poly.standard.polycom_iterative]
    
    References
    ----------
    - [Wikipedia - Horner's method - Efficiency](https://en.wikipedia.org/wiki/Horner%27s_method#Efficiency)
    """
    p = iter(p)
    a0 = next(p, type(x)(0))
    return sumprod_default(p, polyvals(x, start=1), initial=a0, default=MISSING)

def polyval_horner(p, x):
    r"""Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    Uses Horner's method.
    
    `p` must be reversible.
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $\begin{cases}n&n\ge0\\0&n\le0\end{cases}$ scalar additions (`add`) &
    - $\begin{cases}n&n\ge0\\0&n\le0\end{cases}$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`polyval`][poly.standard.polyval]
    - other implementations:
    [`polyval_naive`][poly.standard.polyval_naive],
    [`polyval_iterative`][poly.standard.polyval_iterative]
    - for polynomial arguments: [`polycom_horner`][poly.standard.polycom_horner]
    
    References
    ----------
    - [Wikipedia - Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method)
    """
    p = iter(reversed(p))
    an = next(p, type(x)(0))
    return reduce_default(lambda a, pi: a*x+pi, p, initial=an, default=MISSING)

def polyvals(x, start=0):
    r"""Yield the powers of the value `x`.
    
    $$
        (x^n)_{n\in\mathbb{N}_0} = (1, x, x^2, \dots)
    $$
    
    Uses iterative multiplication to calculate monomials consecutively.
    
    See also
    --------
    - used by: [`polyval_iterative`][poly.standard.polyval_iterative]
    - for polynomial arguments: [`polypows`][poly.standard.polypows]
    """
    if start <= 0:
        yield type(x)(1)
    yield (y := prod_default(repeat(x, start), initial=MISSING, default=x))
    while True:
        y *= x
        yield y

def polyvalzero(p, zero=0):
    """Return the value of polynomial `p` evaluated at point 0.
    
    $$
        p(0)
    $$
    
    More efficient than `polyval(p, 0)`.
    
    Complexity
    ----------
    There are no scalar arithmetic operations.
    
    See also
    --------
    - for any argument: [`polyval`][poly.standard.polyval]
    """
    return next(iter(p), zero)

def polycom(p, q, method='iterative'):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q
    $$
    
    `q` must be a sequence.
    
    Available methods are
    
    - [`naive`][poly.standard.polycom_naive],
    - [`iterative`][poly.standard.polycom_iterative] &
    - [`horner`][poly.standard.polycom_horner] (`p` must be reversible).
    
    See also
    --------
    - implementations: [`polycom_naive`][poly.standard.polycom_naive],
    [`polycom_iterative`][poly.standard.polycom_iterative],
    [`polycom_horner`][poly.standard.polycom_horner]
    - for $q=x-s$: [`polyshift`][poly.standard.polyshift]
    - for scalar arguments: [`polyval`][poly.standard.polyval]
    """
    match method:
        case 'naive':
            return polycom_naive(p, q)
        case 'iterative':
            return polycom_iterative(p, q)
        case 'horner':
            return polycom_horner(p, q)
        case _:
            raise ValueError('Invalid method')

def polycom_naive(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q
    $$
    
    Uses naive repeated multiplication to calculate monomials individually.
    
    `q` must be a sequence.
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be
    
    - $\begin{cases}\frac{n^3m^2}{6}+\frac{n^2m}{2}-\frac{nm^2}{6}-\frac{nm}{2}+n&n\ge0\land m\ge0\\0&n\le0\lor m\le0\end{cases}$ scalar additions (`add`) &
    - $\begin{cases}\frac{n^3m^2}{6}+\frac{n^3m}{6}+n^2m-\frac{nm^2}{6}-\frac{nm}{6}+\frac{n^2}{2}+\frac{n}{2}&n\ge0\land m\ge0\\0&n\le0\lor m\le0\end{cases}$ scalar multiplications (´mul´).
    
    See also
    --------
    - for any implementation: [`polycom`][poly.standard.polycom]
    - other implementations:
    [`polycom_iterative`][poly.standard.polycom_iterative],
    [`polycom_horner`][poly.standard.polycom_horner]
    - for scalar arguments: [`polyval_naive`][poly.standard.polyval_naive]
    """
    p = iter(p)
    a0 = tuple(islice(p, 1)) #(a0,) or polyzero
    return polyadd(a0, *(polyscalarmul(pi, polypow_naive(q, i)) for i, pi in enumerate(p, start=1)))

def polycom_iterative(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q
    $$
    
    Uses iterative multiplication to calculate monomials consecutively.
    
    `q` must be a sequence.
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be
    
    - $\begin{cases}\frac{n^2m^2}{2}+\frac{n^2m}{2}-\frac{nm^2}{2}-\frac{nm}{2}+n&n\ge0\land m\ge0\\0&n\le0\lor m\le0\end{cases}$ scalar additions (`add`) &
    - $\begin{cases}\frac{n^2m^2}{2}+n^2m-\frac{nm^2}{2}+nm+2n-m-1&n\ge1\land m\ge0\\0&n\le1\lor m\le0\end{cases}$ scalar multiplications (´mul´).
    
    See also
    --------
    - for any implementation: [`polycom`][poly.standard.polycom]
    - other implementations:
    [`polycom_naive`][poly.standard.polycom_naive],
    [`polycom_horner`][poly.standard.polycom_horner]
    - uses: [`polypows`][poly.standard.arithmetic.polypows]
    - for scalar arguments: [`polyval_iterative`][poly.standard.polyval_iterative]
    """
    p = iter(p)
    a0 = tuple(islice(p, 1))
    return polyadd(a0, *map(polyscalarmul, p, polypows(q, start=1)))

def polycom_horner(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q
    $$
    
    Uses Horner's method.
    
    `q` must be reversible.
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be
    
    - $\begin{cases}0&m\ge0\lor n<0\\n&m<0\land n>=0\end{cases}$ scalar unary pluses (`pos`),
    - $\begin{cases}\frac{n^2m^2}{2}-\frac{nm^2}{2}+n&n\ge0\land m\ge0\\0&n\le0\lor m\le0\end{cases}$ scalar additions (`add`) &
    - $\begin{cases}\frac{n^2m^2}{2}+\frac{n^2m}{2}-\frac{nm^2}{2}+\frac{nm}{2}+n&n\ge0\land m\ge0\\0&n\le0\lor m\le0\end{cases}$ scalar multiplications (´mul´).
    
    See also
    --------
    - for any implementation: [`polycom`][poly.standard.polycom]
    - other implementations:
    [`polycom_naive`][poly.standard.polycom_naive],
    [`polycom_iterative`][poly.standard.polycom_iterative]
    - for scalar arguments: [`polyval_horner`][poly.standard.polyval_horner]
    """
    p = iter(reversed(p))
    an = tuple(islice(p, 1))
    return reduce_default(lambda a, pi: polyaddc(polymul_naive(a, q), pi), p, initial=an, default=MISSING)

def polyshift(p, s, one=1):
    """Return the polynomial `p` shifted by `s` on the abscissa.
    
    $$
        p(x - s)
    $$
    
    TODO: https://math.stackexchange.com/a/694571/1170417
    
    See also
    --------
    - for polynomial argument: [`polycom`][poly.standard.polycom]
    """
    return polycom_naive(p, (-s, one))

def polyscale(p, a):
    """Return the polynomial `p` scaled by 1/`a` on the abscissa.
    
    $$
        p(ax)
    $$
    
    More efficient than `polycom(p, polyscalarmul(a, polyx))`.
    
    See also
    --------
    - for polynomial argument: [`polycom`][poly.standard.polycom]
    """
    p = iter(p)
    return tuple(chain(islice(p, 1), veclhadamard(p, polyvals(a, start=1))))
