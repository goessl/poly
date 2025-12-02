from itertools import repeat
from vector import vecrshift, vecpos, vecneg, vecadd, vecaddc, vecsub, vecsubc, vecmul, vectruediv, vecfloordiv, vecmod, vecdivmod
from operationcounter import MISSING, reduce_default



__all__ = ('polypos', 'polyneg',
           'polyadd', 'polyaddc', 'polysub', 'polysubc',
           'polyscalarmul', 'polyscalartruediv', 'polyscalarfloordiv',
           'polyscalarmod', 'polyscalardivmod',
           'polymul', 'polymul_naive', 'polymul_karatsuba',
           'polymulx',
           'polypow', 'polypow_naive', 'polypow_binary', 'polypows')



def polypos(p):
    """Return the polynomial with the unary positive operator applied.
    
    $$
        +p
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar unary plus operations (`pos`).
    
    See also
    --------
    - wraps: [`vector.vecpos`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecpos)
    """
    return vecpos(p)

def polyneg(p):
    """Return the polynomial with the unary negative operator applied.
    
    $$
        -p
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar negations (`neg`).
    
    See also
    --------
    - wraps: [`vector.vecneg`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecneg)
    """
    return vecneg(p)

def polyadd(*ps):
    r"""Return the sum of polynomials.
    
    $$
        p_0 + p_1 + \cdots
    $$
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be
    
    - $\min\{n, m\}+1$ scalar additions (`add`).
    
    See also
    --------
    - for constant or monomial summand: [`polyaddc`][poly.standard.polyaddc]
    - wraps: [`vector.vecadd`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecadd)
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polyadd`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polyadd.html)
    """
    return vecadd(*ps)

def polyaddc(p, c, n=0):
    r"""Return the sum of polynomial `p` and a monomial of degree `n`.
    
    $$
        p + cx^n
    $$
    
    More efficient than `polyadd(p, polymonom(n, c))`.
    
    Complexity
    ----------
    There will be
    
    - one scalar addition (`add`) if $n\le\deg p$ or
    - one unary plus operations (`pos`) otherwise.
    
    See also
    --------
    - for polynomial summand: [`polyadd`][poly.standard.polyadd]
    - wraps: [`vector.vecaddc`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecaddc)
    """
    return vecaddc(p, c, i=n)

def polysub(p, q):
    r"""Return the difference of two polynomials.
    
    $$
        p - q
    $$
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be
    
    - $\min\{n, m\}+1$ scalar subtractions (`sub`) &
    - $\begin{cases}m-n&m\ge n\\0&m\le n\end{cases}$ negations (`neg`).
    
    See also
    --------
    - for constant or monomial subtrahend: [`polysubc`][poly.standard.polysubc]
    - wraps: [`vector.vecsub`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecsub)
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polysub`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polysub.html)
    """
    return vecsub(p, q)

def polysubc(p, c, n=0):
    r"""Return the difference of polynomial `p` and a monomial of degree `n`.
    
    $$
        p - cx^n
    $$
    
    More efficient than `polysub(p, polymonom(n, c))`.
    
    Complexity
    ----------
    There will be
    
    - one scalar subtraction (`sub`) if $n\le\deg p$ or
    - one scalar negation (`neg`) otherwise.
    
    See also
    --------
    - for polynomial subtrahend: [`polysub`][poly.standard.polysub]
    - wraps: [`vector.vecsubc`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecsubc)
    """
    return vecsubc(p, c, i=n)

def polyscalarmul(a, p):
    r"""Return the product of a scalar and a polynomial.
    
    $$
        a\cdot p
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar multiplications (`rmul`).
    
    See also
    --------
    - for polynomial factor: [`polymul`][poly.standard.polymul]
    - wraps: [`vector.vecmul`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecmul)
    """
    return vecmul(a, p)

def polyscalartruediv(p, a):
    r"""Return the true division of a polynomial and a scalar.
    
    $$
        \frac{p}{a}
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar true divisions (`truediv`).
    
    See also
    --------
    - wraps: [`vector.vectruediv`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vectruediv)
    """
    return vectruediv(p, a)

def polyscalarfloordiv(p, a):
    r"""Return the floor division of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor x^k
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar floor divisions (`floordiv`).
    
    See also
    --------
    - included in: [`polyscalardivmod`][poly.standard.polyscalardivmod]
    - wraps: [`vector.vecfloordiv`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecfloordiv)
    """
    return vecfloordiv(p, a)

def polyscalarmod(p, a):
    r"""Return the elementwise mod of a polynomial and a scalar.
    
    $$
        \sum_k\left(a_k \mod a \right)x^k
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar modulos (`mod`).
    
    See also
    --------
    - included in: [`polyscalardivmod`][poly.standard.polyscalardivmod]
    - wraps: [`vector.vecmod`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecmod)
    """
    return vecmod(p, a)

def polyscalardivmod(p, a):
    r"""Return the elementwise divmod of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor x^k, \ \sum_k\left(a_k \mod a \right)x^k
    $$
    
    Complexity
    ----------
    For a polynomial of degree $n$ there will be
    
    - $n+1$ scalar divmods (`divmod`).
    
    See also
    --------
    - combines: [`polyscalarfloordiv`][poly.standard.polyscalarfloordiv] & [`polyscalarmod`][poly.standard.polyscalarmod]
    - wraps: [`vector.vecdivmod`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecdivmod)
    """
    return vecdivmod(p, a)

def polymul(*ps, method='naive', one=1):
    r"""Return the product of polynomials.
    
    $$
        p_0 p_1 \cdots
    $$
    
    Available methods are
    
    - [`naive`][poly.standard.polymul_naive] &
    - [`karatsuba`][poly.standard.polymul_karatsuba].
    
    See also
    --------
    - implementations: [`polymul_naive`][poly.standard.arithmetic.polymul_naive],
    [`polymul_karatsuba`][poly.standard.arithmetic.polymul_karatsuba]
    - for scalar factor: [`polyscalarmul`][poly.standard.arithmetic.polyscalarmul]
    - for monomial factor: [`polymulx`][poly.standard.arithmetic.polymulx]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polymul`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polymul.html)
    """
    ps = tuple(map(tuple, ps))
    match method:
        case 'naive':
            return reduce_default(polymul_naive, ps, default=(one,))
        case 'karatsuba':
            return reduce_default(polymul_karatsuba, ps, default=(one,))
        case _:
            raise ValueError('Invalid method')

def polymul_naive(p, q):
    r"""Return the product of two polynomials.
    
    $$
        pq
    $$
    
    Uses naive multiplication and summation.
    
    `q` must be a sequence.
    
    Complexity
    ----------
    For two polynomials of degrees $n$ & $m$ there will be
    
    - $\begin{cases}nm&n\ge1\land m\ge1\\0&n\le0\lor m\le0\end{cases}$ scalar additions (`add`) &
    - $(n+1)(m+1)$ scalar mutliplications (`mul`).
    
    See also
    --------
    - for any implementation: [`polymul`][poly.standard.arithmetic.polymul]
    - other implementations:
    [`polymul_karatsuba`][poly.standard.arithmetic.polymul_karatsuba]
    """
    if not q:
        return () #polyzero
    
    r, sentinel = [], object()
    for i, pi in enumerate(p):
        r.extend([sentinel] * (i+len(q) - len(r)))
        for j, qj in enumerate(q):
            if r[i+j] is sentinel:
                r[i+j] = pi * qj
            else:
                r[i+j] += pi * qj
    return tuple(r)

def _polymul_karatsuba(p, q):
    if len(p) == 1:
        return (p[0] * q[0],)
    
    m = len(p) // 2
    pl, pu = p[:m], p[m:]
    ql, qu = q[:m], q[m:]
    
    z0 = _polymul_karatsuba(pl, ql)
    z2 = _polymul_karatsuba(pu, qu)
    z1 = _polymul_karatsuba(polyadd(pl, pu), polyadd(ql, qu))
    mid = polysub(polysub(z1, z0), z2)
    
    #return vecadd(z0, (0,)*m+mid, (0,)*(2*m)+z2)
    return z0[:m] + polyadd(z0[m:2*m], mid[:m]) + polyadd(z0[2*m:], mid[m:], z2)

def polymul_karatsuba(p, q):
    """Return the product of two polynomials.
    
    $$
        p q
    $$
    
    Uses the Karatsuba algorithm.
    
    Both arguments must be sequences.
    
    TODO: complexity
    
    See also
    --------
    - for any implementation: [`polymul`][poly.standard.arithmetic.polymul]
    - other implementations:
    [`polymul_naive`][poly.standard.arithmetic.polymul_naive]
    
    References
    ----------
    - [Wikipedia - Karatsuba algorithm](https://en.wikipedia.org/wiki/Karatsuba_algorithm)
    """
    if not p or not q:
        return ()
    
    p, q = sorted((q, p), key=len) #p shorter, q longer
    ql, qu = q[:len(p)], q[len(p):]
    rl, ru = _polymul_karatsuba(p, ql), polymul_naive(p, qu)
    #return vecadd(rl, (0,)*len(p)+ru)
    return rl[:len(p)] + polyadd(rl[len(p):], ru)

def polymulx(p, n=1, zero=0):
    """Return the product of polynomial `p` and a monomial of degree `n`.
    
    $$
        px^n
    $$
    
    More efficient than `polymul(p, polymonom(n, c))`.
    
    Complexity
    ----------
    There are no scalar arithmetic operations.
    
    See also
    --------
    - for polynomial factor: [`polymul`][poly.standard.arithmetic.polymul]
    - wraps: [`vector.vecrshift`](https://goessl.github.io/vector/functional/#vector.functional.utility.vecrshift)
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polymulx`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polymulx.html)
    """
    return vecrshift(p, n, zero=zero)

def polypow(p, n, method='naive'):
    """Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    `p` must be a sequence.
    
    Available methods are 
    
    - [`naive`][poly.standard.polypow_naive] &
    - [`binary`][poly.standard.polypow_binary].
    
    TODO: mod parameter
    
    See also
    --------
    - implementations: [`polypow_naive`][poly.standard.arithmetic.polypow_naive],
    [`polypow_binary`][poly.standard.polypow_binary]
    - for sequence of powers: [`polypows`][poly.standard.arithmetic.polypows]
    
    References
    ----------
    - `numpy` equivalent: [`numpy.polynomial.polynomial.polypow`](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.polynomial.polypow.html)
    """
    match method:
        case 'naive':
            return polypow_naive(p, n)
        case 'binary':
            return polypow_binary(p, n)
        case _:
            raise ValueError('Invalid method')

def polypow_naive(p, n, one=1):
    r"""Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Uses repeated multiplication.
    
    `p` must be a sequence.
    
    Complexity
    ----------
    For a polynomial of degree $n$ and exponent $k$ there will be
    
    - $\begin{cases}\frac{n^2k(k-1)}{2}&n\ge0\\0&n\leq0\end{cases}$ scalar additions (`add`) &
    - $\begin{cases}\frac{(nk+2)(n+1)(k-1)}{2}&k>0\\0&k=0\end{cases}$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`polypow`][poly.standard.arithmetic.polypow]
    - other implementations:
    [`polypow_binary`][poly.standard.arithmetic.polypow_binary]
    - uses: [`polypows`][poly.standard.arithmetic.polypows]
    """
    return next(polypows(p, start=n, one=one))

def polypow_binary(p, n, one=1):
    r"""Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Uses exponentiation by squaring.
    
    `p` must be a sequence.
    
    Complexity
    ----------
    For a polynomial of degree $n$ and exponent $k$ let
    
    $$
        \begin{aligned}
            k &= \sum_jb_j2^j \\
            L &= \operatorname{bitlength}(k) \\
            w &= \operatorname{popcount}(k)
        \end{aligned}
    $$
    
    where $j_0$ is the least significant $1$-bit position of $k$,
    $\operatorname{bitlength}(k)$ is the number of bits of the binary
    representation of $k$ and $\operatorname{popcount}(k)$ is the number of
    $1$-bits.
    
    Further define
    
    $$
        \begin{aligned}
            A(k) &= \frac{4^L-1}{3}+\sum_{\substack{i<j \\ b_i=b_j=1}}2^{i+j} \\
            B(k) &= 2(2^L-1)+k-2^{j_0}+\sum_{j\mid b_j=1}2^j\operatorname{popcount}(k\gg(j+1)) \\
            C(k) &= L+w-1.
        \end{aligned}
    $$
    
    ```python
    L = k.bit_length()
    w = k.bit_count()
    C = L + w - 1
    B = 2*(2**L-1) + k-2**((k&-k).bit_length()-1) + sum(2**j * (k>>(j+1)).bit_count() for j in range(L) if ((k>>j)&0x01)==0x01)
    A = (4**L-1)//3 + sum(2**(i+j) for j in range(L) for i in range(j) if ((k>>j)&0x01)==((k>>i)&0x01)==0x01)
    ```
    
    Then there will be
    
    - $A(k)n^2$ scalar additions (`add`) &
    - $A(k)n^2+B(k)n+C(k)$ scalar multiplications (`mul`).
    
    See also
    --------
    - for any implementation: [`polypow`][poly.standard.arithmetic.polypow]
    - other implementations:
    [`polypow_naive`][poly.standard.arithmetic.polypow_naive]
    
    References
    ----------
    - [Wikipedia - Exponentiation by squaring](https://en.wikipedia.org/wiki/Exponentiation_by_squaring)
    - Sequence $C(k)$: [A056791](https://oeis.org/A056791)
    """
    if n == 0:
        return (one,)
    r = None
    while n:
        if n % 2 == 1:
            r = polymul_naive(r, p) if r is not None else p
        p = polymul_naive(p, p)
        n //= 2
    return r

def polypows(p, start=0, one=1):
    r"""Yield the powers of the polynomial `p`.
    
    $$
        (p^n)_{n\in\mathbb{N}_0} = (1, p, p^2, \dots)
    $$
    
    Uses iterative multiplication to calculate powers consecutively.
    
    `p` must be a sequence.
    
    Notes
    -----
    Was first `.evaluation.polycoms` in analogy to
    [`polyvals`][poly.standard.polyvals] but then the submodules `arithmetic` &
    `evaluation` would not be separable (`.evaluation.polycoms` uses
    `.arithmetic.polymul` and `.arithmetic.polypow` uses
    `.evaluation.polycoms`).
    
    See also
    --------
    - used by: [`polypow_naive`][poly.standard.arithmetic.polypow_naive], [`polycom_iterative`][poly.standard.polycom_iterative]
    - for scalar arguments: [`polyvals`][poly.standard.polyvals]
    """
    if start <= 0:
        yield (one,)
    yield (q := reduce_default(polymul_naive, repeat(p, start), initial=MISSING, default=p))
    while True:
        q = polymul_naive(q, p)
        yield q

#def polydiv(n, d):
#    """Return $(q, r)$ such that $n=qd+r$."""
#    #https://en.wikipedia.org/wiki/Polynomial_long_division#Pseudocode
#    n, d, q = vectrim(n), vectrim(d), polyzero
#    while n and polydeg(n)>=polydeg(d):
#        t = vecbasis(len(n)-len(d), n[-1]/d[-1])
#        q = vecadd(q, t)
#        n = vectrim(vecsub(n, polymul(d, t))[:-1])
#    return q, n
