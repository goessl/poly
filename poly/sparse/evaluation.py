from itertools import repeat
from .arithmetic import polysadd, polysscalarmul, polyspow
from operationcounter import sumprod_default



__all__ = ('polysval', 'polysvalzero', 'polyscom', 'polysshift', 'polysscale')



def polysval(p, x):
    """Return the value of polynomial `p` evaluated at point `x`.
    
    $$
        p(x)
    $$
    
    See also
    --------
    - implementations: [`polyval_naive`][poly.standard.polyval_naive],
    [`polyval_iterative`][poly.standard.polyval_iterative],
    [`polyval_horner`][poly.standard.polyval_horner]
    - for consecutive monomials: [`polyvals`][poly.standard.polyvals]
    - for $x=0$: [`polyvalzero`][poly.standard.polyvalzero]
    - for polynomial arguments: [`polycom`][poly.standard.polycom]
    """
    return sumprod_default(p.values(), map(pow, repeat(x), p.keys()), default=type(x)(0))

def polysvalzero(p, zero=0):
    """Return the value of polynomial `p` evaluated at point 0.
    
    $$
        p(0)
    $$
    
    More efficient than `polysval(p, 0)`.
    
    See also
    --------
    - for any argument: [`polysval`][poly.sparse.evaluation.polysval]
    """
    return p.get(0, zero)

def polyscom(p, q):
    r"""Return the polynomial composition of `p` & `q`.
    
    $$
        p\circ q
    $$
    
    See also
    --------
    - for $q=x-s$: [`polysshift`][poly.sparse.evaluation.polysshift]
    - for scalar arguments: [`polysval`][poly.sparse.evaluation.polysval]
    """
    return polysadd(*(polysscalarmul(pi, polyspow(q, i)) for i, pi in p.items()))

def polysshift(p, s, one=1):
    """Return the polynomial `p` shifted by `s` on the abscissa.
    
    $$
        p(x - s)
    $$
    
    See also
    --------
    - for polynomial argument: [`polyscom`][poly.sparse.evaluation.polyscom]
    """
    return polyscom(p, {0:-s, 1:one})

def polysscale(p, a):
    """Return the polynomial `p` scaled by 1/`a` on the abscissa.
    
    $$
        p(ax)
    $$
    
    More efficient than `polyscom(p, polysscalarmul(a, polysx))`.
    
    See also
    --------
    - for polynomial argument: [`polyscom`][poly.sparse.evaluation.polyscom]
    """
    return {i:a**i*pi for i, pi in p.items()}
