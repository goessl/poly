from vector.sparse.utility import vecseq, vecstrim



__all__ = ('polysdeg', 'polyseq', 'polystrim')



def polysdeg(p):
    r"""Return the degree of a polynomial.
    
    $$
        \deg(p)
    $$
    
    Doesn't handle leading zeros, use [`polystrim`][poly.sparse.utility.polystrim]
    if needed.
    
    $\deg(0)=-1$ is used for the empty zero polynomial.
    """
    return max(p.keys(), default=-1)

def polyseq(p, q):
    r"""Return if two polynomials are equal.
    
    $$
        p \overset{?}{=} q
    $$
    
    See also
    --------
    - wraps: [`vecseq`](https://goessl.github.io/vector/sparse/#vector.sparse.utility.vecseq)
    """
    return vecseq(p, q)

def polystrim(p, tol=1e-9):
    """Remove all near zero (`abs(v_i)<=tol`) coefficients."""
    return vecstrim(p, tol=tol)
