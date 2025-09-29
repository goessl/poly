# poly

A Python package for polynomials.
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

This package heavily depends on and is the natural extension of [goessl/vector](https://goessl.github.io/vector).

| Operation               | [poly](functional.md)                                      |
| ----------------------- | ---------------------------------------------------------- |
| **Creation**            |                                                            |
| Zero constant           | [`polyzero`][poly.functional.polyzero]                     |
| One constant            | [`polyone`][poly.functional.polyone]                       |
| Identity                | [`polyx`][poly.functional.polyx]                           |
| Basis                   | [`polymono`][poly.functional.polymono]                     |
| Uniform random          | [`polyrand`][poly.functional.polyrand]                     |
| Normal random           | [`polyrandn`][poly.functional.polyrandn]                   |
| From roots              | [`polyfromroots`][poly.functional.polyfromroots]           |
| **Utility**             |                                                            |
| Comparison              | [`polyeq`][poly.functional.polyeq]                         |
| Trimming                | [`polytrim`][poly.functional.polytrim]                     |
| Rounding                | [`polyround`][poly.functional.polyround]                   |
| Degree                  | [`polydeg`][poly.functional.polydeg]                       |
| **Application**         |                                                            |
| Evaluation              | [`polyval`][poly.functional.polyval]                       |
|                         | [`polyval`][poly.functional.polyval_naive]                 |
|                         | [`polyval`][poly.functional.polyval_iterative]             |
| Evaluation of monomials | [`polyvalgen`][poly.functional.polyvalgen]                 |
| Evaluation at $x=0$     | [`polyvalzero`][poly.functional.polyvalzero]               |
| Composition             | [`polycom`][poly.functional.polycom]                       |
|                         | [`polycom_naive`][poly.functional.polycom_naive]           |
|                         | [`polycom_iterative`][poly.functional.polycom_iterative]   |
|                         | [`polycom_horner`][poly.functional.polycom_horner]         |
| Powers of polynomials   | [`polycomgen`][poly.functional.polycomgen]                 |
| Shift                   | [`polyshift`][poly.functional.polyshift]                   |
| **Arithmetic**          |                                                            |
| Positive                | [`polypos`][poly.functional.polypos]                       |
| Negation                | [`polyneg`][poly.functional.polyneg]                       |
| Addition                | [`polyadd`][poly.functional.polyadd]                       |
| Basis addition          | [`polyaddc`][poly.functional.polyaddc]                     |
| Subtraction             | [`polysub`][poly.functional.polysub]                       |
| Scalar multiplication   | [`polyscalarmul`][poly.functional.polyscalarmul]           |
| Scalar true division    | [`polyscalartruediv`][poly.functional.polyscalartruediv]   |
| Scalar floor division   | [`polyscalarfloordiv`][poly.functional.polyscalarfloordiv] |
| Scalar mod              | [`polyscalarmod`][poly.functional.polyscalarmod]           |
| Scalar divmod           | [`polyscalardivmod`][poly.functional.polyscalardivmod]     |
| Multiplication          | [`polymul`][poly.functional.polymul]                       |
|                         | [`polymul_naive`][poly.functional.polymul_naive]           |
|                         | [`polymul_karatsuba`][poly.functional.polymul_karatsuba]   |
| Multiplication by $x$   | [`polymulx`][poly.functional.polymulx]                     |
| Exponentiation          | [`polypow`][poly.functional.polypow]                       |
|                         | [`polypow_naive`][poly.functional.polypow_naive]           |
|                         | [`polypow_binary`][poly.functional.polypow_binary]         |
| **Calculus**            |                                                            |
| Differentiation         | [`polyder`][poly.functional.polyder]                       |
| Integration             | [`polyantider`][poly.functional.polyantider]               |
| **Conversion**          |                                                            |
| Sympification           | [`polysympify`][poly.functional.polysympify]               |
| Unsympification         | [`polyunsympify`][poly.functional.polyunsympify]           |

## Design

### Coefficient order

Coefficients are stored in ascending order like every sane person would do.

This ways the coefficient indices correspond to the monomial exponent (or the basis index incase of Hermite polynomials).

[The newer `numpy.polynomial` does it this way, the older `numpy.poly1d` does it in reverse](https://numpy.org/doc/stable/reference/routines.polynomials.html#transition-guide). [`sympy` also does it the wrong way](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly.coeffs). 

## Roadmap

- Modules
  - [ ] Hermite polynomials module
  - [ ] Object oriented module
  - [ ] Parallelised module
  - [ ] Multivariate module
- Documentation
  - [ ] [`polymul_karatsuba`][poly.functional.polymul_karatsuba]
  - [ ] Consisten use of `Horner` / `Clenshaw`
- Algorithms
  - [ ] `polysqrt`
  - [ ] `polyroots` only if there is a clean algorithm
  - [ ] [Knuth–Eve evaluation](https://en.wikipedia.org/wiki/Knuth–Eve_algorithm)
  - [ ] [Taylor shift](https://math.stackexchange.com/a/694571/1170417)
  - [ ] Complexity analysis

## License (MIT)

Copyright (c) 2024-2025 Sebastian Gössl

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
