import sympy as sp
from sympy.abc import x



__all__ = ('polysympify', 'polyunsympify')



def polysympify(p, x=x):
    """Return the coefficient iterable `p` as a [`sympy.Poly`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly)."""
    return sp.Poly.from_list(tuple(reversed(tuple(p))), x)

def polyunsympify(p):
    """Return [`sympy.Poly(p)`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly) as a coefficient `tuple`."""
    return tuple(reversed(p.all_coeffs()))
