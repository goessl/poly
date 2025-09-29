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

| Operation               | [poly](functional.md)                                      | [hermite](hermite_functional.md)                                     |
| ----------------------- | ---------------------------------------------------------- | -------------------------------------------------------------------- |
| **Creation**            |                                                            |                                                                      |
| $0$                     | [`polyzero`][poly.functional.polyzero]                     | [`hermzero`][poly.hermite_functional.hermzero]                       |
| $1$                     | [`polyone`][poly.functional.polyone]                       | [`hermone`][poly.hermite_functional.hermone]                         |
| $x$                     | [`polyx`][poly.functional.polyx]                           | [`hermx`][poly.hermite_functional.hermx]                             |
| $x^n$                   | [`polymono`][poly.functional.polymono]                     | [`hermmono`][poly.hermite_functional.hermmono]                       |
| $H_0$                   |                                                            | [`H0`][poly.hermite_functional.H0]                                   |
| $H_1$                   |                                                            | [`H1`][poly.hermite_functional.H1]                                   |
| $H_2$                   |                                                            | [`H2`][poly.hermite_functional.H2]                                   |
| $H_n$                   |                                                            | [`herm`][poly.hermite_functional.herm]                               |
| Uniform random          | [`polyrand`][poly.functional.polyrand]                     | [`hermrand`][poly.hermite_functional.hermrandn]                      |
| Normal random           | [`polyrandn`][poly.functional.polyrandn]                   | [`hermrandn`][poly.hermite_functional.hermrandn]                     |
| From roots              | [`polyfromroots`][poly.functional.polyfromroots]           | [`hermfromroots`][poly.hermite_functional.hermfromroots]             |
| **Utility**             |                                                            |                                                                      |
| Comparison              | [`polyeq`][poly.functional.polyeq]                         | [`hermeq`][poly.hermite_functional.hermeq]                           |
| Trimming                | [`polytrim`][poly.functional.polytrim]                     | [`hermtrim`][poly.hermite_functional.hermtrim]                       |
| Rounding                | [`polyround`][poly.functional.polyround]                   | [`hermround`][poly.hermite_functional.hermround]                     |
| Degree                  | [`polydeg`][poly.functional.polydeg]                       | [`hermdeg`][poly.hermite_functional.hermdeg]                         |
| **Conversion**          |                                                            |                                                                      |
| $H\to P$                |                                                            | [`herm2poly`][poly.hermite_functional.herm2poly]                     |
|                         |                                                            | [`herm2poly`][poly.hermite_functional.herm2poly_naive]               |
|                         |                                                            | [`herm2poly`][poly.hermite_functional.herm2poly_iterative]           |
|                         |                                                            | [`herm2poly`][poly.hermite_functional.herm2poly_clenshaw]            |
| $P\to H$                |                                                            | [`poly2herm`][poly.hermite_functional.poly2herm]                     |
|                         |                                                            | [`poly2herm_naive`][poly.hermite_functional.poly2herm_naive]         |
|                         |                                                            | [`poly2herm_iterative`][poly.hermite_functional.poly2herm_iterative] |
|                         |                                                            | [`poly2herm_horner`][poly.hermite_functional.poly2herm_horner]       |
| **Evaluation**          |                                                            |                                                                      |
| Evaluation              | [`polyval`][poly.functional.polyval]                       | [`hermval`][poly.hermite_functional.hermval]                         |
|                         | [`polyval_naive`][poly.functional.polyval_naive]           | [`hermval_naive`][poly.hermite_functional.hermval_naive]             |
|                         | [`polyval_iterative`][poly.functional.polyval_iterative]   | [`hermval_iterative`][poly.hermite_functional.hermval_iterative]     |
| Evaluation of basis     | [`polyvalgen`][poly.functional.polyvalgen]                 | [`hermvalgen`][poly.hermite_functional.hermvalgen]                   |
| Evaluation at $x=0$     | [`polyvalzero`][poly.functional.polyvalzero]               | [`hermvalzero`][poly.hermite_functional.hermvalzero]                 |
| Composition             | [`polycom`][poly.functional.polycom]                       |                                                                      |
|                         | [`polycom_naive`][poly.functional.polycom_naive]           |                                                                      |
|                         | [`polycom_iterative`][poly.functional.polycom_iterative]   |                                                                      |
|                         | [`polycom_horner`][poly.functional.polycom_horner]         |                                                                      |
| Powers of polynomials   | [`polycomgen`][poly.functional.polycomgen]                 |                                                                      |
| Shift                   | [`polyshift`][poly.functional.polyshift]                   |                                                                      |
| **Arithmetic**          |                                                            |                                                                      |
| Positive                | [`polypos`][poly.functional.polypos]                       | [`hermpos`][poly.hermite_functional.hermpos]                         |
| Negation                | [`polyneg`][poly.functional.polyneg]                       | [`hermneg`][poly.hermite_functional.hermneg]                         |
| Addition                | [`polyadd`][poly.functional.polyadd]                       | [`hermadd`][poly.hermite_functional.hermadd]                         |
| Basis addition          | [`polyaddc`][poly.functional.polyaddc]                     | [`hermaddc`][poly.hermite_functional.hermaddc]                       |
| Subtraction             | [`polysub`][poly.functional.polysub]                       | [`hermsub`][poly.hermite_functional.hermsub]                         |
| Scalar multiplication   | [`polyscalarmul`][poly.functional.polyscalarmul]           | [`hermscalarmul`][poly.hermite_functional.hermscalarmul]             |
| Scalar true division    | [`polyscalartruediv`][poly.functional.polyscalartruediv]   | [`hermscalartruediv`][poly.hermite_functional.hermscalartruediv]     |
| Scalar floor division   | [`polyscalarfloordiv`][poly.functional.polyscalarfloordiv] | [`hermscalarfloordiv`][poly.hermite_functional.hermscalarfloordiv]   |
| Scalar mod              | [`polyscalarmod`][poly.functional.polyscalarmod]           | [`hermscalarmod`][poly.hermite_functional.hermscalarmod]             |
| Scalar divmod           | [`polyscalardivmod`][poly.functional.polyscalardivmod]     | [`hermscalardivmod`][poly.hermite_functional.hermscalardivmod]       |
| Multiplication          | [`polymul`][poly.functional.polymul]                       | [`hermmul`][poly.hermite_functional.hermmul]                         |
|                         | [`polymul_naive`][poly.functional.polymul_naive]           | [`hermmul_naive`][poly.hermite_functional.hermmul_naive]             |
|                         | [`polymul_karatsuba`][poly.functional.polymul_karatsuba]   |                                                                      |
| Multiplication by $x$   | [`polymulx`][poly.functional.polymulx]                     | [`hermmulx`][poly.hermite_functional.hermmulx]                       |
| Exponentiation          | [`polypow`][poly.functional.polypow]                       | [`hermpow`][poly.hermite_functional.hermpow]                         |
|                         | [`polypow_naive`][poly.functional.polypow_naive]           | [`hermpow_naive`][poly.hermite_functional.hermpow_naive]             |
|                         | [`polypow_binary`][poly.functional.polypow_binary]         | [`hermpow_binary`][poly.hermite_functional.hermpow_binary]           |
| **Calculus**            |                                                            |                                                                      |
| Differentiation         | [`polyder`][poly.functional.polyder]                       | [`hermder`][poly.hermite_functional.hermder]                         |
| Integration             | [`polyantider`][poly.functional.polyantider]               | [`hermantider`][poly.hermite_functional.hermantider]                 |
| **Conversion**          |                                                            |                                                                      |
| Sympification           | [`polysympify`][poly.functional.polysympify]               | [`hermsympify`][poly.hermite_functional.hermsympify]                 |
| Unsympification         | [`polyunsympify`][poly.functional.polyunsympify]           |                                                                      |

## Design

### Coefficient order

Coefficients are stored in ascending order like every sane person would do.

This ways the coefficient indices correspond to the monomial exponent (or the basis index incase of Hermite polynomials).

[The newer `numpy.polynomial` does it this way, the older `numpy.poly1d` does it in reverse](https://numpy.org/doc/stable/reference/routines.polynomials.html#transition-guide). [`sympy` also does it the wrong way](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly.coeffs). 

## Roadmap

- Modules
  - [x] Hermite polynomials module
  - [ ] Object oriented module
  - [ ] Parallelised module
  - [ ] Multivariate module
- Documentation
  - [ ] [`polymul_karatsuba`][poly.functional.polymul_karatsuba]
  - [ ] Consistent `See also` links; graph view?
  - [ ] Consistent use of `Horner` / `Clenshaw`
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
