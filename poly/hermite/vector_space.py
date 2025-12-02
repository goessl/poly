from vector import vecpos, vecneg, vecadd, vecaddc, vecsub, vecsubc, vecmul, vectruediv, vecfloordiv, vecmod, vecdivmod



__all__ = ('hermpos', 'hermneg', 'hermadd', 'hermaddc', 'hermsub', 'hermsubc',
           'hermscalarmul', 'hermscalartruediv', 'hermscalarfloordiv',
           'hermscalarmod', 'hermscalardivmod')



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
    - in standard monomial basis: [`polypos`][poly.standard.vector_space.polypos]
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
    - in standard monomial basis: [`polyneg`][poly.standard.vector_space.polyneg]
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
    - for single coefficient summand ($cH_n$): [`hermaddc`][poly.hermite.vector_space.hermaddc]
    - in standard monomial basis: [`polyadd`][poly.standard.vector_space.polyadd]
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
    - for series summand: [`hermadd`][poly.hermite.vector_space.hermadd]
    - in standard monomial basis: [`polyaddc`][poly.standard.vector_space.polyaddc]
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
    - for single coefficient subtrahend ($cH_n$): [`hermaddc`][poly.hermite.vector_space.hermaddc]
    - in standard monomial basis: [`polysub`][poly.standard.vector_space.polysub]
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
    - for series minuend: [`hermsub`][poly.hermite.vector_space.hermsub]
    - in standard monomial basis: [`polysubc`][poly.standard.vector_space.polysubc]
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
    - in standard monomial basis: [`polyscalarmul`][poly.standard.vector_space.polyscalarmul]
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
    - in standard monomial basis: [`polyscalartruediv`][poly.standard.vector_space.polyscalartruediv]
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
    - in standard monomial basis: [`polyscalarfloordiv`][poly.standard.vector_space.polyscalarfloordiv]
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
    - in standard monomial basis: [`polyscalarmod`][poly.standard.vector_space.polyscalarmod]
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
    - in standard monomial basis: [`polyscalardivmod`][poly.standard.vector_space.polyscalardivmod]
    - wraps: [`vector.vecdivmod`](https://goessl.github.io/vector/functional/#vector.functional.vector_space.vecdivmod)
    """
    return vecdivmod(h, a)
