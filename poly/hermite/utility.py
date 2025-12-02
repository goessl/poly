from vector import veceq, vectrim, vecround



__all__ = ('hermeq', 'hermtrim', 'hermround', 'hermdeg')



def hermeq(g, h):
    r"""Return if two Hermite polynomial series are equal.
    
    $$
        g \overset{?}{=} h
    $$
    
    See also
    --------
    - standard monomial basis: [`polyeq`][poly.standard.polyeq]
    - wraps: [`vector.veceq`](https://goessl.github.io/vector/functional/#vector.functional.veceq)
    """
    return veceq(g, h)

def hermtrim(h, tol=1e-9):
    r"""Remove all leading near zero (`abs(h_i)<=tol`) coefficients.
    
    $$
        \sum_{k=0}^nh_kH_k \ \text{where} \ n=\max\{\, k\mid |a_k|>\text{tol}\,\}\cup\{-1\}
    $$
    
    See also
    --------
    - standard monomial basis: [`polytrim`][poly.standard.polytrim]
    - wraps: [`vector.vectrim`](https://goessl.github.io/vector/functional/#vector.functional.vectrim)
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermtrim`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermtrim.html)
    """
    return vectrim(h, tol)

def hermround(h, ndigits=None):
    r"""Round all coefficients to the given precision.
    
    $$
        \sum_k\text{round}_\text{ndigits}(h_k)H_k
    $$
    
    See also
    --------
    - standard monomial basis: [`polyround`][poly.standard.polyround]
    - wraps: [`vector.vecround`](https://goessl.github.io/vector/functional/#vector.functional.vecround)
    """
    return vecround(h, ndigits)

def hermdeg(h):
    r"""Return the degree of a Hermite polynomial series.
    
    $$
        \deg(h)
    $$
    
    $\deg(0)=-1$ is used for the empty zero Hermite polynomial series.
    
    Doesn't handle leading zeros, use
    [`hermtrim`][poly.hermite.hermtrim] if needed.
    
    Notes
    -----
    $\deg(0)=-\infty$ is more commonly used but the expected return type is an
    `int` and `-math.inf` is of type `float`. Therefore the $\deg(0)=-1$
    convention was choosen to keep the return type consistent.
    
    See also
    --------
    - standard monomial basis: [`polydeg`][poly.standard.polydeg]
    
    References
    ----------
    - [Wikipedia - Degree of a polynomial - Degree of the zero polynomial](https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial)
    """
    return sum(1 for _ in h) - 1
