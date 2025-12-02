from math import comb
from itertools import count, islice, repeat, tee
from .hilbert_space import hermweightis
from vector import vecbasis, vecrshift, vecpos, vecneg, vecadd, vecladd, vecaddc, vecsub, vecsubc, vecmul, vectruediv, vecltruediv, vecfloordiv, vecmod, vecdivmod, veclhadamard
from operationcounter import MISSING, reduce_default



__all__ = ('hermpos', 'hermneg', 'hermadd', 'hermaddc', 'hermsub', 'hermsubc',
           'hermscalarmul', 'hermscalartruediv', 'hermscalarfloordiv',
           'hermscalarmod', 'hermscalardivmod',
           'hermmul', 'hermmul_naive', 'hermmulx', 'hermmulHn',
           'hermpow', 'hermpow_naive', 'hermpow_binary', 'hermmulpow', 'hermpows')



def hermpos(h):
    """Return the Hermite polynomial series with the unary positive operator applied.
    
    $$
        +h
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar unary plus operations (`pos`).
    
    See also
    --------
    - wraps: [`vector.vecpos`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecpos)
    """
    return vecpos(h)

def hermneg(h):
    """Return the Hermite polynomial series with the unary negative operator applied.
    
    $$
        -h
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar negations (`neg`).
    
    See also
    --------
    - wraps: [`vector.vecneg`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecneg)
    """
    return vecneg(h)

def hermadd(*hs):
    r"""Return the sum of Hermite polynomial series.
    
    $$
        h_0 + h_1 + \cdots
    $$
    
    Complexity
    ----------
    For two Hermite polynomial series of degrees $n$ & $m$ there will be
    
    - $\min\{n, m\}+1$ scalar additions (`add`).
    
    See also
    --------
    - for single coefficient summand ($cH_n$): [`hermaddc`][poly.hermite.arithmetic.hermaddc]
    - wraps: [`vector.vecadd`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecadd)
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermadd`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermadd.html)
    """
    return vecadd(*hs)

def hermaddc(h, c, n=0):
    r"""Return the sum of Hermite polynomial series `h` and the `n`-th Hermite polynomial.
    
    $$
        h + cH_n
    $$
    
    More efficient than `hermadd(h, vecbasis(n, c))`.
    
    Complexity
    ----------
    There will be
    
    - one scalar addition (`add`) if $n\le\deg h$ or
    - one unary plus operations (`pos`) otherwise.
    
    See also
    --------
    - for series summand: [`hermadd`][poly.hermite.arithmetic.hermadd]
    - wraps: [`vector.vecaddc`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecaddc)
    """
    return vecaddc(h, c, n)

def hermsub(g, h):
    r"""Return the difference of two Hermite polynomial series.
    
    $$
        g - h
    $$
    
    Complexity
    ----------
    For two Hermite polynomial series of degrees $n$ & $m$ there will be
    
    - $\min\{n, m\}+1$ scalar subtractions (`sub`) &
    - $\begin{cases}m-n&m\ge n\\0&m\le n\end{cases}$ negations (`neg`).
    
    See also
    --------
    - for single coefficient subtrahend ($cH_n$): [`hermaddc`][poly.hermite.arithmetic.hermaddc]
    - wraps: [`vector.vecsub`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecsub)
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermsub`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermsub.html)
    """
    return vecsub(g, h)

def hermsubc(h, c, n=0):
    r"""Return the difference of Hermite polynomial series `h` and the `n`-th Hermite polynomial.
    
    $$
        h - cH_n
    $$
    
    More efficient than `hermsub(h, vecbasis(n, c))`.
    
    Complexity
    ----------
    There will be
    
    - one scalar subtraction (`sub`) if $n\le\deg h$ or
    - one scalar negation (`neg`) otherwise.
    
    See also
    --------
    - for series minuend: [`hermsub`][poly.hermite.arithmetic.hermsub]
    - wraps: [`vector.vecsubc`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecsubc)
    """
    return vecsubc(h, c, n)

def hermscalarmul(a, h):
    r"""Return the product of a scalar and a Hermite polynomial series.
    
    $$
        a\cdot h
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar multiplications (`rmul`).
    
    See also
    --------
    - wraps: [`vector.vecmul`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecmul)
    """
    return vecmul(a, h)

def hermscalartruediv(h, a):
    r"""Return the true division of a Hermite polynomial series and a scalar.
    
    $$
        \frac{h}{a}
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar true divisions (`truediv`).
    
    See also
    --------
    - wraps: [`vector.vectruediv`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vectruediv)
    """
    return vectruediv(h, a)

def hermscalarfloordiv(h, a):
    r"""Return the floor division of a Hermite polynomial series and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor H_k
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar floor divisions (`floordiv`).
    
    See also
    --------
    - wraps: [`vector.vecmfloordiv`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecmfloordiv)
    """
    return vecfloordiv(h, a)

def hermscalarmod(h, a):
    r"""Return the elementwise mod of a Hermite polynomial series and a scalar.
    
    $$
        \sum_k\left(a_k \mod a \right)H_k
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar modulos (`mod`).
    
    See also
    --------
    - wraps: [`vector.vecmod`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecmod)
    """
    return vecmod(h, a)

def hermscalardivmod(h, a):
    r"""Return the elementwise divmod of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor H_k, \ \sum_k\left(a_k \mod a \right)H_k
    $$
    
    Complexity
    ----------
    For a Hermite polynomial series of degree $n$ there will be
    
    - $n+1$ scalar divmods (`divmod`).
    
    See also
    --------
    - wraps: [`vector.vecdivmod`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecdivmod)
    """
    return vecdivmod(h, a)

def hermmul(*hs, method='naive', one=1):
    r"""Return the product of Hermite polynomial series.
    
    $$
        h_0 h_1 \cdots
    $$
    
    Available methods are
    
    - [`naive`][poly.hermite.hermmul_naive].
    
    TODO: Karatsuba possible?
    
    See also
    --------
    - implementations: [`hermmul_naive`][poly.hermite.arithmetic.hermmul_naive]
    - for monomial factor: [`hermmulx`][poly.hermite.arithmetic.hermmulx]
    - for Hermite polynomial factor: [`hermmulHn`][poly.hermite.arithmetic.hermmulHn]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermmul`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermmul.html)
    """
    hs = tuple(map(tuple, hs))
    match method:
        case 'naive':
            return reduce_default(hermmul_naive, hs, default=(one,))
        case _:
            raise ValueError('Invalid method')

def hermmul_naive(g, h):
    r"""Return the product of two Hermite polynomials series.
    
    $$
        g h
    $$
    
    Uses naive multiplication and summation.
    
    `h` must be a sequence.
    
    Complexity
    ----------
    For two Hermite polynomial series of degrees $n$ & $m$ (let $N=\min\{n, m\}$, $M=\max\{n, m\}$) there will be
    
    - $\begin{cases}\frac{(N+1)(N+2)(3M+3-N)}{6}-(n+m+1) & n\geq0\land m\geq0 \\ 0 & n<0\lor m<0 \end{cases}$ scalar additions (`add`) &
    - $\begin{cases}\frac{(N+1)(N+2)(3M+3-N)}{3} & n\geq0\land m\geq0 \\ 0 & n<0\lor m<0 \end{cases}$ scalar multiplications (`mul`).
    
    Notes
    -----
    From [1] eq. A.8 we know that
    
    $$
        \tilde{H}_i\tilde{H}_j = \sum_{k=0}^{\min\{i,j\}}k!\binom{i}{k}\binom{j}{k}\tilde{H}_{i+j-2k}.
    $$
    
    Therefore we have
    
    $$
        \begin{aligned}
            &H_i(x)H_j(x) &&\mid H_j(x) = 2^\frac{j}{2}\tilde{H}_j(\sqrt{2}x) \\
            &= 2^\frac{i+j}{2} \tilde{H}_i(\sqrt{2}x) \tilde{H}_j(\sqrt{2}x) &&\mid \tilde{H}_i\tilde{H}_j = \sum_{k=0}^{\min\{i, j\}}k!\binom{i}{k}\binom{j}{k}\tilde{H}_{i+j-2k} \\
            &= 2^\frac{i+j}{2} \sum_{k=0}^{\min\{i, j\}} k!\binom{i}{k}\binom{j}{k} \tilde{H}_{i+j-2k}(\sqrt{2}x) &&\mid \cdot1=2^{k-k} \\
            &= \sum_{k=0}^{\min\{i, j\}} 2^kk!\binom{i}{k}\binom{j}{k}2^\frac{i+j-2k}{2} \tilde{H}_{i+j-2k}(\sqrt{2}x) &&\mid H_j(x)=2^\frac{j}{2} \tilde{H}_j(\sqrt{2}x) \\
            &= \sum_{k=0}^{\min\{i, j\}}2^kk!\binom{i}{k}\binom{j}{k} H_{i+j-2k}(x)
        \end{aligned}
    $$
    
    See also
    --------
    - for any implementation: [`hermmul`][poly.hermite.arithmetic.hermmul]
    
    References
    ----------
    - [1] Zhi-yuan Huang & Jia-an Yan: Introduction to Infinite Dimensional Stochastic Analysis. [10.1007/978-94-011-4108-6](https://doi.org/10.1007/978-94-011-4108-6)
    - [Wikipedia - Hermite polynomials - Definition](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition)
    """
    if not h:
        return () #hermzero
    
    r, sentinel = [], object()
    for i, gi in enumerate(g):
        r.extend([sentinel] * (i+len(h) - len(r)))
        for j, hj in enumerate(h):
            for k, f in enumerate(islice(hermweightis(), min(i, j)+1)):
                if r[i+j-2*k] is sentinel:
                    r[i+j-2*k] = f * comb(i, k)*comb(j, k) * gi*hj
                else:
                    r[i+j-2*k] += f * comb(i, k)*comb(j, k) * gi*hj
    return tuple(r)

def hermmulx(h, zero=0):
    r"""Return the product of Hermite polynomial series `h` and a monomial of degree `n`.
    
    $$
        h(x)x^n
    $$
    
    More efficient than `hermmul(h, hermx)`.
    
    Complexity
    ----------
    For a Hermite function series of degree $n$ there will be
    
    - $\begin{cases}n-1 & n<0 \\ 0 & n\leq0 \end{cases}$ scalar additions (`add`),
    - $\begin{cases}n-1 & n<0 \\ 0 & n\leq0 \end{cases}$ scalar multiplications (`mul`) &
    - $n+1$ scalar divisions (`truediv`) by two.
    
    Notes
    --------
    For a single Hermite polynomial
    
    $$
        \begin{aligned}
            H_{k+1}(x) &= 2xH_k(x)-H_k'(x) \\
            H_k'(x) &= 2kH_{k-1}(x) \\
            \Rightarrow \quad xH_k(x) &= \frac{H_{k+1}(x)}{2}+kH_{k-1}(x)
        \end{aligned}
    $$
    
    and therefore for a series
    
    $$
        \begin{aligned}
            xh(x) &= x\sum_kh_kH_k(x) \\
            &= \sum_kh_kxH_k(x) \\
            &= \sum_kh_k\left(\frac{1}{2}H_{k+1}(x)+kH_{k-1}(x)\right) \\
            &= \sum_k\frac{h_k}{2}H_{k+1}(x)+\sum_kkh_kH_{k-1}(x) \\
            &= \sum_k\frac{h_{k-1}}{2}H_k(x)+\sum_k(k+1)h_{k+1}H_k(x) \\
            &= \sum_k\left(\frac{h_{k-1}}{2}+(k+1)h_{k+1}\right)H_k(x)
        \end{aligned}
    $$
    
    See also
    --------
    - for Hermite polynomial series as factor: [`hermmul`][poly.hermite.arithmetic.hermmul]
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermmulx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermmulx.html)
    
    References
    ----------
    - [Wikipedia - Hermite polynomials - Recurrence relation](https://en.wikipedia.org/wiki/Hermite_polynomials#Recurrence_relation)
    """
    im1, ip1 = tee(h)
    _, h1 = next(ip1, None), next(ip1, zero)
    return vecrshift(vecladd(vecltruediv(im1, 2), veclhadamard(count(2), ip1)), 1, zero=h1)

def hermmulHn(h, n, zero=0):
    r"""Return the product of Hermite polynomial series `h` and the `n`-th Hermite polynomial.
    
    $$
        hH_n
    $$
    
    More efficient than `hermmul(h, vecbasis(n))`.
    
    TODO: Complexity
    
    See also
    --------
    - for a Hermite polynomial series factor: [`hermmul`][poly.hermite.arithmetic.hermmul]
    """
    r = []
    for i, hi in enumerate(h):
        r.extend([zero] * (i+n - len(r)+1))
        for k, f in enumerate(islice(hermweightis(), min(i, n)+1)):
            if r[i+n-2*k] is zero:
                r[i+n-2*k] = f * comb(i, k)*comb(n, k) * hi
            else:
                r[i+n-2*k] += f * comb(i, k)*comb(n, k) * hi
    return tuple(r)

def hermpow(h, n, method='binary'):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    `h` must be a sequence.
    
    Available methods are 
    
    - [`naive`][poly.hermite.hermpow_naive] &
    - [`binary`][poly.hermite.hermpow_binary].
    
    See also
    --------
    - implementations: [`hermpow_naive`][poly.hermite.hermpow_naive],
    [`hermpow_binary`][poly.hermite.hermpow_binary],
    - `numpy` equivalent: [`numpy.polynomial.hermite.hermpow`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.hermite.hermpow.html)
    """
    match method:
        case 'naive':
            return hermpow_naive(h, n)
        case 'binary':
            return hermpow_binary(h, n)
        case _:
            raise ValueError('Invalid method')

def hermpow_naive(h, n, one=1):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    Uses repeated multiplication.
    
    `h` must be a sequence.
    
    See also
    --------
    - for any implementation: [`hermpow`][poly.hermite.hermpow]
    - other implementations: [`hermpow_binary`][poly.hermite.hermpow_binary]
    """
    return next(hermpows(h, start=n, one=one))

def hermpow_binary(h, n, one=1):
    """Return the Hermite polynomial series `h` raised to the nonnegative `n`-th power.
    
    $$
        h^n
    $$
    
    Uses exponentiation by squaring.
    
    `h` must be a sequence.
    
    See also
    --------
    - for any implementations: [`hermpow`][poly.hermite.hermpow]
    - other implementations[`hermpow_naive`][poly.hermite.hermpow_naive]
    
    References
    ----------
    - [Wikipedia - Exponentiation by squaring](https://en.wikipedia.org/wiki/Exponentiation_by_squaring)
    """
    if n==0:
        return (one,) #hermone
    r = None
    while n:
        if n % 2 == 1:
            r = hermmul(r, h) if r is not None else h
        h = hermmul(h, h)
        n //= 2
    return r

def hermmulpow(alpha, one=1):
    r"""Return product of Hermite Polynomials.
    
    $$
        H^\alpha = \prod_iH_i^{\alpha_i}
    $$
    """
    r = (one,) #hermone
    for i, ai in enumerate(alpha):
        r = hermmul(r, hermpow(vecbasis(i), ai, one=one))
    return r

def hermpows(h, start=0, one=1):
    r"""Yield the powers of the Hermite polynomial series `h`.
    
    $$
        (h^n)_{n\in\mathbb{N}_0} = (1, h, h^2, \dots)
    $$
    
    Uses iterative multiplication to calculate powers consecutively.
    
    `h` must be a sequence.
    
    See also
    --------
    - used by: [`hermpow_naive`][poly.hermite.arithmetic.hermpow_naive]
    - for scalar arguments: [`hermvals`][poly.hermite.hermvals]
    """
    if start <= 0:
        yield (one,) #hermone
    yield (g := reduce_default(hermmul_naive, repeat(h, start), initial=MISSING, default=h))
    while True:
        g = hermmul_naive(g, h)
        yield g
