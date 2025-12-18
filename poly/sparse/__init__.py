"""Sparse polynomials.

```python
>>> from poly import polysval
>>> p = {0:1, 3:2} #p(x)=2*x^3+1
>>> polysval(p, 5)

```

Prefixed by `polys...` (poly - sparse).

All functions accept polynomials and return them as `dict`s (monomial-exponent:coefficient).

Index keys are expected to be integers.

## Docstring conventions

Summary

Math notation (vector notation if possible, index notation, domain & codomain)

More information ("More efficient than ...").

Complexity
----------
For a polynomial with $n$ elements there will be - $x$ scalar additions (`add`), ...

Notes
-----
Design choices

See also
--------
Similar functions

References
----------
Wikipedia, numpy, ...
"""

from .creation import *
from .conversion import *
from .utility import *
from .evaluation import *
from .arithmetic import *
from .calculus import *
