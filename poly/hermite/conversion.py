from itertools import count
from ..standard import polyzero, polyadd, polyaddc, polysub, polyscalarmul, polymulx
from .arithmetic import hermadd, hermaddc, hermscalarmul, hermmulx
from .creation import hermzero, herms, herm_explicit, hermmono, hermmonos
from vector import vecdot
import sympy as sp
from sympy.abc import x



__all__ = ('herm2poly', 'herm2poly_naive', 'herm2poly_iterative', 'herm2poly_clenshaw',
           'poly2herm', 'poly2herm_naive', 'poly2herm_iterative', 'poly2herm_horner',
           'hermsympify')



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
    return polyadd(*map(polyscalarmul, h, herms()))

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
    return hermadd(*map(hermscalarmul, p, hermmonos()))

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

def hermsympify(h, x=x):
    """Return the Hermite polynomial series `h` as a [`sympy.Poly`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly)."""
    return vecdot(h, (sp.hermite_poly(i, x, polys=True) for i in count()), zero=sp.Poly(0, x))
