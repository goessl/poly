from itertools import count
from vector import vecdot
import sympy as sp
from sympy.abc import x



__all__ = ('hermsympify', )



def hermsympify(h, x=x):
    """Return the Hermite polynomial series `h` as a [`sympy.Poly`](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly)."""
    return vecdot(h, (sp.hermite_poly(i, x, polys=True) for i in count()), zero=sp.Poly(0, x))
