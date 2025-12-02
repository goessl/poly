from math import sqrt, pi, factorial
from itertools import count
from vector import vecdot, vecabsq



__all__ = ('hermweight', 'hermweights', 'hermweighti', 'hermweightis',
           'hermabs', 'hermabsq', 'hermabsqi',
           'hermdot', 'hermdoti')



def hermweight(n):
    r"""Return the normalisation factor of Hermite polynomials.
    
    $$
        \int_\mathbb{R}H_n^2(x)e^{-x^2}\,\mathrm{d}x = \sqrt{\pi}2^nn!
    $$
    
    See also
    --------
    - for the integral factor: [`hermweighti`][poly.hermite.hermweighti]
    - for all normalisation factors: [`hermweights`][poly.hermite.hermweights]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    return sqrt(pi) * hermweighti(n)

def hermweights(start=0):
    r"""Yield the normalisation factors of Hermite polynomials.
    
    $$
        \begin{gathered}
            \left(\int_\mathbb{R}H_n^2(x)e^{-x^2}\,\mathrm{d}x\right)_{n\in\mathbb{N}_0} = \sqrt{\pi}(2^nn!)_{n\in\mathbb{N}_0} \\
            = \sqrt{\pi}(1, 2, 8, 48, 384, 3840, 46080, 645120, \dots)
        \end{gathered}
    $$
    
    See also
    --------
    - for the integral factors: [`hermweightis`][poly.hermite.hermweightis]
    - for a single normalisation factor: [`hermweight`][poly.hermite.hermweight]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    sqrtpi = sqrt(pi)
    for f in hermweightis(start):
        yield sqrtpi * f

def hermweighti(n):
    r"""Return the integral factor of the normalisation factor of Hermite polynomials.
    
    $$
        \frac{1}{\sqrt{\pi}}\int_\mathbb{R}H_n^2(x)e^{-x^2}\,\mathrm{d}x = 2^nn!
    $$
    
    See also
    --------
    - for the full normalisation factor: [`hermweight`][poly.hermite.hermweight]
    - for all normalisation factors: [`hermweightis`][poly.hermite.hermweightis]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    return 2**n * factorial(n)

def hermweightis(start=0):
    r"""Yield the integral factors of the normalisation factors of Hermite polynomials.
    
    $$
        \begin{gathered}
            \frac{1}{\sqrt{\pi}}\left(\int_\mathbb{R}H_n^2(x)e^{-x^2}\,\mathrm{d}x\right)_{n\in\mathbb{N}_0} = (2^nn!)_{n\in\mathbb{N}_0} \\
            = (1, 2, 8, 48, 384, 3840, 46080, 645120, \dots)
        \end{gathered}
    $$
    
    See also
    --------
    - for all full normalisation factors: [`hermweights`][poly.hermite.hermweights]
    - for a single normalisation factor: [`hermweighti`][poly.hermite.hermweighti]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    - [A000165](https://oeis.org/A000165)
    """
    yield (f := 2**start * factorial(start))
    for n in count(start+1):
        f *= 2 * n
        yield f

def hermabs(h, conjugate=False, zero=0):
    r"""Return the norm of a Hermite polynomial series.
    
    $$
        ||h||_H = \sqrt{\int_\mathbb{R}h^{(*)}(x)h(x)e^{-x^2}\,\mathrm{d}x} = \sqrt{\sqrt{\pi}\sum_k2^kk!h^{(*)}_kh_k}
    $$
    
    See also
    --------
    - for the norm squared: [`hermabsq`][poly.hermite.hermabsq]
    - similar to: [`vector.vecabs`](https://goessl.github.io/vector/functional/#vector.functional.hilbert_space.vecabs)
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    return hermabsq(h, conjugate=conjugate, zero=zero)**0.5

def hermabsq(h, conjugate=False, zero=0):
    r"""Return squared norm of a Hermite polynomial series.
    
    $$
        ||h||_H^2 = \int_\mathbb{R}h^{(*)}(x)h(x)e^{-x^2}\,\mathrm{d}x = \sqrt{\pi}\sum_k2^kk!h^{(*)}_kh_k
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar conjugations (`conjugate`) (if selected),
    - $n+1$/$2(n+1)$ scalar multiplications (`mul`) without/with weights &
    - $\begin{cases}n&n\ge1\\0&n\le1\end{cases}$ scalar additions (`add`).
    
    See also
    --------
    - for the norm: [`hermabs`][poly.hermite.hermabs]
    - for the integral factor: [`hermabsqi`][poly.hermite.hermabsqi]
    - similar to: [`vector.vecabsq`](https://goessl.github.io/vector/functional/#vector.functional.hilbert_space.vecabsq)
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    return sqrt(pi) * hermabsqi(h, conjugate=conjugate, zero=zero)

def hermabsqi(h, conjugate=False, zero=0):
    r"""Return the integral factor of the norm squared of a Hermite polynomial series.
    
    $$
        \frac{1}{\sqrt{\pi}}||h||_H^2 = \frac{1}{\sqrt{\pi}}\int_\mathbb{R}h^{(*)}(x)h(x)e^{-x^2}\,\mathrm{d}x = \sum_k2^kk!h^{(*)}_kh_k
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar conjugations (`conjugate`) (if selected),
    - $\begin{cases}n&n\ge0\\0&n\le0\end{cases}$ scalar additions (`add`) &
    - $2(n+1)$ scalar multiplications (`mul`) without/with weights.
    
    See also
    --------
    - for the full norm squared: [`hermabsq`][poly.hermite.hermabsq]
    - wraps: [`vector.vecabsq`](https://goessl.github.io/vector/functional/#vector.functional.hilbert_space.vecabsq)
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    return vecabsq(h, weights=hermweightis(), conjugate=conjugate, zero=zero)

def hermdot(g, h, conjugate=False, zero=0):
    r"""Return the inner product with respect to the Hermite polynomial weight function.
    
    $$
        \left<g\mid h\right>_H = \int_\mathbb{R}g^{(*)}(x)h(x)e^{-x^2}\,\mathrm{d}x = \sqrt{\pi}\sum_k2^kk!g^{(*)}_kh_k
    $$
    
    Complexity
    ----------
    For two Hermite polynomial series of degrees $n$ & $m$ there will be
    
    - $\min\{n, m\}+1$ scalar conjugations (`conjugate`) (if selected),
    - $\begin{cases}n&n\ge0\\0&n\le0\end{cases}$ scalar additions (`add`) &
    - $2(n+1)$ scalar multiplications (`mul`) without/with weights.
    
    See also
    --------
    - for the integral factor: [`hermdoti`][poly.hermite.hermdoti]
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Orthogonality](https://en.wikipedia.org/wiki/Hermite_polynomials#Orthogonality)
    """
    return sqrt(pi) * hermdoti(g, h, conjugate=conjugate, zero=zero)

def hermdoti(g, h, conjugate=False, zero=0):
    r"""Return the integral factor of the inner product with respect to the Hermite polynomial weight function.
    
    $$
        \frac{1}{\sqrt{\pi}}\left<g\mid h\right>_H = \frac{1}{\sqrt{\pi}}\int_\mathbb{R}g^{(*)}(x)h(x)e^{-x^2}\,\mathrm{d}x = \sum_k2^kk!g^{(*)}_kh_k
    $$
    
    See also
    --------
    - for the full inner product: [`hermdot`][poly.hermite.hermdot]
    - wraps: [`vector.vecdot`](https://goessl.github.io/vector/functional/#vector.functional.hilbert_space.vecdot)
    """
    return vecdot(g, h, weights=hermweightis(), conjugate=conjugate, zero=zero)
