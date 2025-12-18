from math import perm
from .evaluation import polysval
from .arithmetic import polysaddc



__all__ = ('polysder', 'polysantider')



def polysder(p, k=1):
    """Return the `k`-th derivative of polynomial `p`.
    
    $$
        p^{(k)}
    $$
    """
    return {i-k:perm(i, k)*pi for i, pi in p.items() if i-k>=0}

def polysantider(p, c=0, b=0):
    r"""Return the antiderivative of polynomial `p`.
    
    $$
        \int_b^xp(y)\,\mathrm{d}y+c
    $$
    """
    P = {i+1:pi/(i+1) for i, pi in p.items()}
    return polysaddc(P, c-polysval(P, b) if b else c)
