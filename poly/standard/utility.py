from vector import veclen, veceq, vectrim, vecround



__all__ = ('polydeg', 'polyeq', 'polytrim', 'polyround')



def polydeg(p):
    r"""Return the degree of a polynomial.
    
    $$
        \deg(p)
    $$
    
    Doesn't handle leading zeros, use [`polytrim`][poly.standard.polytrim]
    if needed.
    
    $\deg(0)=-1$ is used for the empty zero polynomial.
    
    Notes
    -----
    $\deg(0)=-\infty$ is more commonly used but the expected return type is an
    `int` where `-math.inf` is of type `float`. Therefore the $\deg(0)=-1$
    convention was choosen to keep the return type consistent.
    
    See also
    --------
    - wraps: [`vector.veclen`](https://goessl.github.io/vector/functional/#vector.functional.utility.veclen)
    
    References
    ----------
    - [Wikipedia - Degree of a polynomial - Degree of the zero polynomial](https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial)
    """
    return veclen(p) - 1

def polyeq(p, q):
    r"""Return if two polynomials are equal.
    
    $$
        p \overset{?}{=} q
    $$
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be at most
    
    - $\min\{n, m\}+1$ scalar comparisons (`eq`) &
    - $|n-m|$ scalar boolean evaluations (`bool`).
    
    See also
    --------
    - wraps: [`vector.veceq`](https://goessl.github.io/vector/functional/#vector.functional.utility.veceq)
    """
    return veceq(p, q)

def polytrim(p, tol=1e-9):
    r"""Remove all leading near zero (`abs(a_i)<=tol`) coefficients.
    
    $$
        \sum_{k=0}^na_kx^k \ \text{where} \ n=\max\{\, k\mid |a_k|>\text{tol}\,\}\cup\{-1\}
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar absolute evaluations (`abs`) &
    - $n+1$ scalar comparisons (`gt`).
    
    See also
    --------
    - wraps: [`vector.vectrim`](https://goessl.github.io/vector/functional/#vector.functional.utility.vectrim)
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polytrim`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polytrim.html)
    """
    return vectrim(p, tol=tol)

def polyround(p, ndigits=None):
    r"""Round all coefficients to the given precision.
    
    $$
        \sum_k\text{round}_\text{ndigits}(a_k)x^k
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar roundings (`round`).
    
    See also
    --------
    - wraps: [`vector.vecround`](https://goessl.github.io/vector/functional/#vector.functional.utility.vecround)
    """
    return vecround(p, ndigits)
