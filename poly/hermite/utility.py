from vector import veclen, veceq, vectrim



__all__ = ('hermdeg', 'hermeq', 'hermtrim')



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
    - wraps: [`vector.veclen`](https://goessl.github.io/vector/functional/#vector.functional.utility.veclen)
    
    References
    ----------
    - [Wikipedia - Degree of a polynomial - Degree of the zero polynomial](https://en.wikipedia.org/wiki/Degree_of_a_polynomial#Degree_of_the_zero_polynomial)
    """
    return veclen(h) - 1

def hermeq(g, h):
    r"""Return if two Hermite polynomial series are equal.
    
    $$
        g \overset{?}{=} h
    $$
    
    See also
    --------
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
    - wraps: [`vector.vectrim`](https://goessl.github.io/vector/functional/#vector.functional.vectrim)
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermtrim`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermtrim.html)
    """
    return vectrim(h, tol)
