from math import perm
from itertools import count, islice
from .evaluation import hermval_iterative
from .arithmetic import hermaddc
from vector import vecrshift, vechadamard, veclhadamardtruediv



__all__ = ('hermder', 'hermantider')



def hermder(h, k=1):
    r"""Return the `k`-th derivative of Hermite polynomial series `h`.
    
    $$
        h^{(k)}
    $$
    
    Notes
    -----
    For Hermite polynomials:
    
    $$
        \begin{aligned}
            \frac{d}{dx}H_n(x) &= 2nH_{n-1}(x) \\
            \frac{d^2}{dx^2}H_n(x) &= 2^2n(n-1)H_{n-2}(x) \\
            &\vdots \\
            \frac{d^k}{dx^k}H_n(x) &= 2^kn(n-1)\cdots(n-k+1)H_{n-k}(x) \\
            &= 2^k(n)_kH_{n-k}(x) \\
            &= 2^k\frac{n!}{(n-k)!}H_{n-k}(x) \\
            &= 2^k{}_nP_kH_{n-k}(x)
        \end{aligned}
    $$
    
    And for Hermite polynomial series:
    
    $$
        \begin{aligned}
            \frac{d^k}{dx^k}h(x) &= \frac{d^k}{dx^k}\sum_{l=0}^nh_lH_l(x) \\
            &= \sum_{l=k}^nh_l2^k{}_lP_kH_{l-k}(x) \\
            &= 2^k\sum_{l=0}^{n-k}h_{l+k}\,{}_{l+k}P_kH_l(x)
        \end{aligned}
    $$
    
    Where $(n)_k$ is the [Falling factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials#Properties)
    and ${}_nP_k$ is the [number of k-permutations of n](https://en.wikipedia.org/wiki/Permutation#k-permutations_of_n)
    with $(n)_k=\frac{n!}{(n-k)!}={}_nP_k$ (falling factorials are
    used because their definition appears in the derivation; permutations are
    used because a fast implementation is provided by `math.perm`).
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    
    See also
    --------
    - in standard monomial basis: [`polyder`][poly.standard.polyder]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermder`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermder.html)
    """
    tk = 2**k
    return vechadamard((tk*perm(l+k, k) for l in count()), islice(h, k, None))

def hermantider(h, b=0, c=0):
    r"""Return the antiderivative of Hermite polynomial series `h`.
    
    $$
        \int_b^xh(y)\,\mathrm{d}y + c
    $$
    
    TODO: higher antiderivatives, complexity
    
    Notes
    -----
    Let
    
    $$
        I_{b,c}[f](x) = \int_b^xf(y)\,\mathrm{d}y+c.
    $$
    
    Then we have for Hermite polynomials:
    
    $$
        \begin{aligned}
            I_{b,c}[H_n](x) &= \int_b^xH_n(y)\,\mathrm{d}y+c \\
            &\qquad \mid H_n'=2nH_{n-1} \quad \Leftrightarrow \quad H_n^{(-1)}=\frac{H_{n+1}}{2(n+1)}+C \\
            &= \left.\frac{H_{n+1}}{2(n+1)}\right|_b^x+c \\
            &= \frac{H_{n+1}(x)}{2(n+1)}-\frac{H_{n+1}(b)}{2(n+1)}+c
        \end{aligned}
    $$
    
    For Hermite polynomials series:
    
    $$
        \begin{aligned}
            I_{b,c}[h](x) &= I_{b,c}\left[\sum_kh_kH_k\right](x) \\
            &= \int_b^x\sum_kh_kH_k(y)\,\mathrm{d}y+c \\
            &= \sum_kh_k\int_b^xH_k(y)\,\mathrm{d}y+c \\
            &= \sum_kh_k\left(\frac{H_{k+1}(x)}{2(k+1)}-\frac{H_{k+1}(b)}{2(k+1)}\right)+c \\
            &= \sum_kh_{k-1}\left(\frac{H_k(x)}{2k}-\frac{H_k(b)}{2k}\right)+c \\
            &= \sum_k\frac{h_{k-1}}{2k}H_k(x)-\sum_k\frac{h_{k-1}}{2k}H_k(b)+c
        \end{aligned}
    $$
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    
    See also
    --------
    - in standard monomial basis: [`polyantider`][poly.standard.polyantider]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermint`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermint.html)
    """
    H = vecrshift(veclhadamardtruediv(h, count(2, 2)), 1)
    return hermaddc(H, c-hermval_iterative(H, b))

r'''
def hermintegral(h, a=1):
    r"""Return the integral factor of the integral over a Hermite polynomial series.
    
    $$
        \sqrt{\frac{a}{\pi}}\int_\mathbb{R}e^{-ax^2}h(x)\,\mathrm{d}x
    $$
    
    Notes
    -----
    $$
        \begi{aligned}
            \int_\mathbb{R}e^{-ax^2}h(x)\,\mathrm{d}x \\
            &= \int_\mathbb{R}e^{-\left(\sqrt{a}x\right)^2}h(x)\,\mathrm{d}x \\
        \end{aligned}
    $$
    """
    return sum((h2n* factorial(2*n) * (Fraction(1, a)-1)**n / factorial(n) for n, h2n in enumerate(islice(h, 0, None, 2))), start=Fraction(0))

def hermpolyintegral(p, a=1):
    r"""Return $\sqrt{\frac{a}{\pi}}\int_\mathbb{R}e^{-ax^2}p(x)dx$."""
    return sum(((Fraction(factorial2(2*n-1, exact=True), (2*a)**n) * p2n if n>0 else p2n) for n, p2n in enumerate(islice(p, 0, None, 2))), start=Fraction(0))
'''
