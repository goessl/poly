# poly

A native polynomial Python module.
```python
>>> from poly import polyval
>>> from fractions import Fraction
>>> p = (Fraction(1, 2), Fraction(3, 4))
>>> x = Fraction(5, 6)
>>> polyval(p, x)
Fraction(9, 8)
```

## Installation

```
pip install git+https://github.com/goessl/poly.git
```

## Usage

This package provides functions, similar to the ones found in `numpy.polynomial.polynomial`, but in pure Python.
This means it can handle Pythons native unlimited `int`s, `Fraction`s and every object that implements the needed scalar arithmetic without `numpy`s lossy conversion.

To represent polynomials sequences like `list`s and `tuple`s are used (coefficients in ascending order). The results are always returned as `tuple`s. The arguments are consumed but not altered.
The functions signatures try to mimic `numpy.polynomial.polynomial`, most often with less parameters because of less functionality.
Operations that are associative, like `polyadd` and `polymul`, allow arbitrary many arguments, even none.

If actual `polyzero`s (`()`, not some untrimmed polynomial like `(0,)`) are passed as arguments, the functions also correctly return actual `polyzero`s (`()`) or `polyone`s (`(1,)`).
The user is responsible to `polytrim` arguments and results.

Provided are fundamental polynomials:
 - `polyzero`: $0$,
 - `polyone`: $1$,
 - `polyx`: $x$,
 - `polymono(n, c=1)`: $c\cdot x^n$;

evaluation:
 - `polyval(p, x)`: $p(x)$,
 - `polycom(p, q)`: $p(q)$;

arithmetic operations:
 - `polyadd(*ps)`: $p_0 + p_1 + \dots$,
 - `polysub(p, q)`: $p - q$,
 - `polymul(*ps)`: $p_0 \cdot p_1 \cdot \dots$,
 - `polydiv(n, d)`: `(q, r)` such that $n=qd+r$,
 - `polypow(p, n)`: $p^n$;

calculus:
 - `polyder(p)`: $p'$,
 - `polyint(p, c=0)`: $\int p(x)dx$

## License (MIT)

Copyright (c) 2024 Sebastian Gössl

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
