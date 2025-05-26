from vector import Vector
from ._polyfunctions import polyzero, polyval, polyder, polyint, polysympify



__all__ = ['Poly']



class Poly(Vector):
    @property
    def deg(self):
        """Return the degree of this polynomial."""
        return len(self) - 1
    
    def __call__(self, x):
        return polyval(self, x)
    
    
    def __mul__(self, other):
        if isinstance(other, Poly):
            return Poly(polymul(self, other))
        return vecmul(self, other)
    
    def __pow__(self, other):
        return Poly(polypow(self, other))
    
    
    def der(self, n=1):
        return Poly(polyder(self, n))
    
    def antider(self, n=1):
        return Poly(polyint(self, n))
    
    
    #python stuff
    def _sympy_(self):
        return polysympify(self)
    
    def __str__(self):
        return ''.join(f'{i:+}x^{i}' for i, hi in enumerate(self))
    
    def _repr_latex_(self):
        return '$' + ''.join(f'{i:+}x^{{{i}}}' for i, hi in enumerate(self)) + '$'



Poly.ZERO = Poly(polyzero)
"""Zero polynomial."""
