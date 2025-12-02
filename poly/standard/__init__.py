r"""Polynomials in standard monomial basis.

$$
    p(x)=\sum_{k=0}^na_kx^k \qquad P=(x_n)_{n\in\mathbb{N}_0}=(1, x, x^2, \dots)
$$

Prefixed by `poly...` (polynomial).

All functions accept **single exhaustible iterables**, if not stated otherwise.

If the result of a function is a polynomial, it is returned as a `tuple`.

Associative operations, like [`polyadd`][poly.standard.polyadd] and
[`polymul`][poly.standard.polymul], allow arbitrary many
arguments, even none.

The functions are **type-independent**. However, the data types used must
*support necessary scalar operations*. For instance, for polynomial addition,
components must be addable.

## Docstring convention

Summary

Math notation

Complexity
----------
For a polynomial of degree $n$ there will be

- $x$ scalar additions (`add`),
- $y$ scalar subtractions, ...

Order: `pos`, `neg`, `add`, `sub`, `mul`

Notes
-----
design choices

See also
--------
other implementations, delegations/wrappings, ...

References
----------

`numpy` equivalents, Wikipedia, ...)
"""



from .creation import *
from .utility import *
from .evaluation import *
from .arithmetic import *
from .calculus import *
from .conversion import *
