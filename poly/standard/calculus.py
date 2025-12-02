from math import perm
from itertools import count, islice
from .evaluation import polyval_iterative
from vector import vecrshift, vechadamard, vechadamardtruediv



__all__ = ('polyder', 'polyantider')



def polyder(p, k=1):
    r"""Return the `k`-th derivative of polynomial `p`.
    
    $$
        p^{(k)}
    $$
    
    Complexity
    ----------
    For the $k$-th derivative of a polynomial of degree $n$ there will be
    
    - $\begin{cases}n-k+1&k\le n\\0&k>n\end{cases}$ scalar multiplications with integers (`rmul`).
    
    Notes
    -----
    For monomials:
    
    $$
        \begin{aligned}
            \frac{d}{dx}x^n &= nx^{n-1} \\
            \frac{d^2}{dx^2}x^n &= n(n-1)x^{n-2} \\
            &\vdots \\
            \frac{d^k}{dx^k}x^n &= n(n-1)\cdots(n-k+1)x^{n-k} \\
            &= (n)_kx^{n-k} \\
            &= \frac{n!}{(n-k)!}x^{n-k} \\
            &= {}_nP_kx^{n-k}
        \end{aligned}
    $$
    
    And for polynomials:
    
    $$
        \begin{aligned}
            \frac{d^k}{dx^k}p(x) &= \frac{d^k}{dx^k}\sum_{l=0}^np_lx^l \\
            &= \sum_{l=k}^np_l\,{}_lP_kx^{l-k} \\
            &= \sum_{l=0}^{n-k}p_{l+k}\,{}_{l+k}P_kx^l
        \end{aligned}
    $$
    
    Where $(n)_k$ is the [Falling factorial](https://en.wikipedia.org/wiki/Falling_and_rising_factorials#Properties)
    and ${}_nP_k$ is the [number of k-permutations of n](https://en.wikipedia.org/wiki/Permutation#k-permutations_of_n)
    with $(n)_k=\frac{n!}{(n-k)!}={}_nP_k$ (falling factorials are
    used because their definition appears in the derivation; permutations are
    used because a fast implementation is provided by `math.perm`).
    
    References
    ----------
    - [more_itertools.polynomial_derivative](https://more-itertools.readthedocs.io/en/stable/api.html#more_itertools.polynomial_derivative)
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polyder`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyder.html)
    - [Wikipedia - Potenzregel - Höhere Ableitung einer Potenzfunktion mit natürlichem Exponenten](https://de.wikipedia.org/wiki/Potenzregel#H%C3%B6here_Ableitung_einer_Potenzfunktion_mit_nat%C3%BCrlichem_Exponenten)
    """
    return vechadamard((perm(l+k, k) for l in count()), islice(p, k, None))

def polyantider(p, c=0, b=0):
    r"""Return the antiderivative of polynomial `p`.
    
    $$
        \int_b^xp(y)\,\mathrm{d}y+c
    $$
    
    TODO: Higher antiderivatives, complexity
    
    Notes
    -----
    Let
    
    $$
        I_{b, c}[f](x) = \int_b^xf(y)\,\mathrm{d}y+c.
    $$
    
    Then we have for monomials:
    
    $$
        \begin{aligned}
            I_{b,c}[\cdot^n](x) &= \int_b^xy^n\,\mathrm{d}y+c \\
            &= \left.\frac{y^{n+1}}{n+1}\right|_{y=b}^x+c \\
            &= \frac{x^{n+1}}{n+1}-\frac{b^{n+1}}{n+1}+c
        \end{aligned}
    $$
    
    For polynomials:
    
    $$
        \begin{aligned}
            I_{b,c}[p](x) &= I_{b,c}\left[\sum_ka_k\cdot^k\right](x) \\
            &= \int_b^x\sum_ka_ky^k\,\mathrm{d}y+c \\
            &= \sum_ka_k\left(\frac{x^{k+1}}{k+1}-\frac{b^{k+1}}{k+1}\right)+c \\
            &= \sum_ka_{k-1}\left(\frac{x^k}{k}-\frac{b^k}{k}\right)+c \\
            &= \sum_k\frac{a_{k-1}}{k}x^k-\sum_k\frac{a_{k-1}}{k}b^k+c
        \end{aligned}
    $$
    
    Notes
    -----
    Integration is called antiderivative (`antider`)
    instead of integrate (`int`) to avoid keyword collisions.
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polyint`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyint.html)
    """
    P = vechadamardtruediv(p, count(1))
    return vecrshift(P, 1, zero=c-b*polyval_iterative(P, b) if b else c)
