from itertools import repeat
from vector.sparse import vecsrshift, vecspos, vecsneg, vecsadd, vecsaddc, vecssub, vecssubc, vecsmul, vecstruediv, vecsfloordiv, vecsmod, vecsdivmod
from operationcounter import MISSING, reduce_default



__all__ = ('polyspos', 'polysneg', 'polysadd', 'polysaddc', 'polyssub', 'polyssubc',
           'polysscalarmul', 'polysscalartruediv', 'polysscalarfloordiv', 'polysscalarmod', 'polysscalardivmod',
           'polysmul', 'polysmul_naive', 'polysmulx', 'polyspow', 'polyspow_naive', 'polyspow_binary', 'polyspows')



def polyspos(p):
    """Return the polynomial with the unary positive operator applied.
    
    $$
        +p
    $$
    
    See also
    --------
    - wraps: [`vecspos`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecspos)
    """
    return vecspos(p)

def polysneg(p):
    """Return the polynomial with the unary negative operator applied.
    
    $$
        -p
    $$
    
    See also
    --------
    - wraps: [`vecsneg`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsneg)
    """
    return vecsneg(p)

def polysadd(*ps):
    r"""Return the sum of polynomials.
    
    $$
        p_0+p_1+\cdots
    $$
    
    See also
    --------
    - wraps: [`vecsadd`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsadd)
    """
    return vecsadd(*ps)

def polysaddc(p, c, n=0):
    """Return the sum of polynomial `p` and a monomial of degree `n`.
    
    $$
        p+cx^n
    $$
    
    See also
    --------
    - wraps: [`vecsaddc`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsaddc)
    """
    return vecsaddc(p, c=c, i=n)

def polyssub(p, q):
    r"""Return the difference of two polynomials.
    
    $$
        p-q
    $$
    
    See also
    --------
    - wraps: [`vecssub`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecssub)
    """
    return vecssub(p, q)

def polyssubc(p, c, n=0):
    """Return the difference of polynomial `p` and a monomial of degree `n`.
    
    $$
        p-cx^n
    $$
    
    See also
    --------
    - wraps: [`vecssubc`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecssubc)
    """
    return vecssubc(p, c=c, i=n)

def polysscalarmul(a, p):
    """Return the product of a scalar and a polynomial.
    
    $$
        ap
    $$
    
    See also
    --------
    - wraps: [`vecsmul`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsmul)
    """
    return vecsmul(a, p)

def polysscalartruediv(p, a):
    r"""Return the true division of a polynomial and a scalar.
    
    $$
        \frac{p}{a}
    $$
    
    See also
    --------
    - wraps: [`vecstruediv`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecstruediv)
    """
    return vecstruediv(p, a)

def polysscalarfloordiv(p, a):
    r"""Return the floor division of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor x^k
    $$
    
    See also
    --------
    - wraps: [`vecsfloordiv`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsfloordiv)
    """
    return vecsfloordiv(p, a)

def polysscalarmod(p, a):
    r"""Return the elementwise mod of a polynomial and a scalar.
    
    $$
        \sum_k\left(a_k \mod a\right)x^k
    $$
    
    See also
    --------
    - wraps: [`vecsmod`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsmod)
    """
    return vecsmod(p, a)

def polysscalardivmod(p, a):
    r"""Return the elementwise divmod of a polynomial and a scalar.
    
    $$
        \sum_k\left\lfloor\frac{a_k}{a}\right\rfloor x^k, \ \sum_k\left(a_k \mod a\right)x^k
    $$
    
    See also
    --------
    - wraps: [`vecsdivmod`](https://goessl.github.io/vector/sparse/#vector.sparse.vector_space.vecsdivmod)
    """
    return vecsdivmod(p, a)

def polysmul(*ps, method='naive', one=1):
    r"""Return the product of polynomials.
    
    $$
        p_0 p_1 \cdots
    $$
    
    Available methods are
    
    - [`naive`][poly.sparse.arithmetic.polysmul_naive].
    
    See also
    --------
    - implementations: [`polysmul_naive`][poly.sparse.arithmetic.polysmul_naive]
    - for scalar factor: [`polysscalarmul`][poly.sparse.arithmetic.polysscalarmul]
    - for monomial factor: [`polysmulx`][poly.sparse.arithmetic.polysmulx]
    """
    match method:
        case 'naive':
            return reduce_default(polysmul_naive, ps, default={0:one})
        case _:
            raise ValueError('Invalid method')

def polysmul_naive(p, q):
    """Return the product of two polynomials.
    
    $$
        pq
    $$
    
    Uses naive multiplication and summation.
    
    See also
    --------
    - for any implementation: [`polysmul`][poly.sparse.arithmetic.polysmul]
    """
    r = {}
    for i, pi in p.items():
        for j, qj in q.items():
            if i+j not in r:
                r[i+j] = pi * qj
            else:
                r[i+j] += pi * qj
    return r

def polysmulx(p, n=1):
    """Return the product of polynomial `p` and a monomial of degree `n`.
    
    $$
        px^n
    $$
    
    More efficient than `polysmul(p, polysmonom(n))`.
    
    Complexity
    ----------
    There are no scalar arithmetic operations.
    
    See also
    --------
    - for polynomial factor: [`polysmul`][poly.sparse.arithmetic.polysmul]
    - wraps: [`vector.vecsrshift`](https://goessl.github.io/vector/sparse/#vector.sparse.utility.vecsrshift)
    """
    return vecsrshift(p, n)

def polyspow(p, n, method='naive'):
    """Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Available methods are 
    
    - [`naive`][poly.sparse.arithmetic.polyspow_naive] &
    - [`binary`][poly.sparse.arithmetic.polyspow_binary].
    
    TODO: mod parameter
    
    See also
    --------
    - implementations: [`polyspow_naive`][poly.sparse.arithmetic.polyspow_naive],
    [`polyspow_binary`][poly.sparse.arithmetic.polyspow_binary]
    - for sequence of powers: [`polyspows`][poly.sparse.arithmetic.polyspows]
    """
    match method:
        case 'naive':
            return polyspow_naive(p, n)
        case 'binary':
            return polyspow_binary(p, n)
        case _:
            raise ValueError('Invalid method')

def polyspow_naive(p, n, one=1):
    r"""Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Uses repeated multiplication.
    
    See also
    --------
    - for any implementation: [`polyspow`][poly.sparse.arithmetic.polyspow]
    - other implementations:
    [`polyspow_binary`][poly.sparse.arithmetic.polyspow_binary]
    - uses: [`polyspows`][poly.sparse.arithmetic.polyspows]
    """
    return next(polyspows(p, start=n, one=one))

def polyspow_binary(p, n, one=1):
    r"""Return the polynomial `p` raised to the nonnegative `n`-th power.
    
    $$
        p^n
    $$
    
    Uses exponentiation by squaring.
    
    See also
    --------
    - for any implementation: [`polyspow`][poly.sparse.arithmetic.polyspow]
    - other implementations:
    [`polyspow_naive`][poly.sparse.arithmetic.polyspow_naive]
    
    References
    ----------
    - [Wikipedia - Exponentiation by squaring](https://en.wikipedia.org/wiki/Exponentiation_by_squaring)
    - Sequence $C(k)$: [A056791](https://oeis.org/A056791)
    """
    if n == 0:
        return {0:one} #polysone
    r = None
    while n:
        if n % 2 == 1:
            r = polysmul_naive(r, p) if r is not None else p
        p = polysmul_naive(p, p)
        n //= 2
    return r

def polyspows(p, start=0, one=1):
    r"""Yield the powers of the polynomial `p`.
    
    $$
        (p^n)_{n\in\mathbb{N}_0} = (1, p, p^2, \dots)
    $$
    
    Uses iterative multiplication to calculate powers consecutively.
    
    See also
    --------
    - used by: [`polyspow_naive`][poly.sparse.arithmetic.polyspow_naive], [`polyscom_iterative`][poly.sparse.evaluation.polyscom]
    - for scalar arguments: [`polyvals`][poly.standard.evaluation.polyvals]
    """
    if start <= 0:
        yield {0:one} #polysone
    yield (q := reduce_default(polysmul_naive, repeat(p, start), initial=MISSING, default=p))
    while True:
        q = polysmul_naive(q, p)
        yield q
