from vector.sparse.conversion import vecstod, vecdtos
from sympy import Poly
from sympy.abc import x



__all__ = ('polystod', 'polydtos', 'polyssympify', 'polysunsympify')



def polystod(p, zero=0):
    """Return a sparse polynomial (`dict`) as a dense polynomial (`tuple`)."""
    return vecstod(p, zero=zero)

def polydtos(p):
    """Return a dense polynomial (`tuple`) as a sparse polynomial (`dict`).
    
    The resulting `dict` is not [trimmed][poly.sparse.utility.polystrim].
    """
    return vecdtos(p)

def polyssympify(p, gen=x):
    """Return the coefficient mapping `p` as a [`sympy.Poly`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly)."""
    return Poly.from_dict(p, gen)

def polysunsympify(p, gen=x):
    """Return [`sympy.Poly(p)`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly) as a coefficient `dict`."""
    return {i[0]:pi for i, pi in p.as_dict(native=True).items()}
