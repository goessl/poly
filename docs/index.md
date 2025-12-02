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

| Operation                  | [poly](standard.md)                                                 | [hermite](hermite.md)                                     |
| -------------------------- | ------------------------------------------------------------------- | --------------------------------------------------------- |
| **Creation**               |                                                                     |                                                           |
| $0$                        | [`polyzero`][poly.standard.creation.polyzero]                       | [`hermzero`][poly.hermite.hermzero]                       |
| $1$                        | [`polyone`][poly.standard.creation.polyone]                         | [`hermone`][poly.hermite.hermone]                         |
| $x$                        | [`polyx`][poly.standard.creation.polyx]                             | [`hermx`][poly.hermite.hermx]                             |
| $x^n$                      | [`polymono`][poly.standard.creation.polymono]                       | [`hermmono`][poly.hermite.hermmono]                       |
| $(x^n)_{n\in\mathbb{N}_0}$ | [`polymonos`][poly.standard.creation.polymono]                      | [`hermmono`][poly.hermite.hermmono]                       |
| $H_0$                      |                                                                     | [`H0`][poly.hermite.H0]                                   |
| $H_1$                      |                                                                     | [`H1`][poly.hermite.H1]                                   |
| $H_2$                      |                                                                     | [`H2`][poly.hermite.H2]                                   |
| $H_n$                      |                                                                     | [`herm`][poly.hermite.herm]                               |
| Uniform random             | [`polyrand`][poly.standard.creation.polyrand]                       | [`hermrand`][poly.hermite.hermrandn]                      |
| Normal random              | [`polyrandn`][poly.standard.creation.polyrandn]                     | [`hermrandn`][poly.hermite.hermrandn]                     |
| From roots                 | [`polyfromroots`][poly.standard.creation.polyfromroots]             | [`hermfromroots`][poly.hermite.hermfromroots]             |
| **Utility**                |                                                                     |                                                           |
| Degree                     | [`polydeg`][poly.standard.utility.polydeg]                          | [`hermdeg`][poly.hermite.hermdeg]                         |
| Comparison                 | [`polyeq`][poly.standard.utility.polyeq]                            | [`hermeq`][poly.hermite.hermeq]                           |
| Trimming                   | [`polytrim`][poly.standard.utility.polytrim]                        | [`hermtrim`][poly.hermite.hermtrim]                       |
| Rounding                   | [`polyround`][poly.standard.utility.polyround]                      | [`hermround`][poly.hermite.hermround]                     |
| **Conversion**             |                                                                     |                                                           |
| $H\to P$                   |                                                                     | [`herm2poly`][poly.hermite.herm2poly]                     |
|                            |                                                                     | [`herm2poly`][poly.hermite.herm2poly_naive]               |
|                            |                                                                     | [`herm2poly`][poly.hermite.herm2poly_iterative]           |
|                            |                                                                     | [`herm2poly`][poly.hermite.herm2poly_clenshaw]            |
| $P\to H$                   |                                                                     | [`poly2herm`][poly.hermite.poly2herm]                     |
|                            |                                                                     | [`poly2herm_naive`][poly.hermite.poly2herm_naive]         |
|                            |                                                                     | [`poly2herm_iterative`][poly.hermite.poly2herm_iterative] |
|                            |                                                                     | [`poly2herm_horner`][poly.hermite.poly2herm_horner]       |
| **Evaluation**             |                                                                     |                                                           |
| Evaluation                 | [`polyval`][poly.standard.evaluation.polyval]                       | [`hermval`][poly.hermite.hermval]                         |
|                            | [`polyval_naive`][poly.standard.evaluation.polyval_naive]           | [`hermval_naive`][poly.hermite.hermval_naive]             |
|                            | [`polyval_iterative`][poly.standard.evaluation.polyval_iterative]   | [`hermval_iterative`][poly.hermite.hermval_iterative]     |
| Evaluation of basis        | [`polyvals`][poly.standard.evaluation.polyvals]                     | [`hermvals`][poly.hermite.hermvals]                       |
| Evaluation at $x=0$        | [`polyvalzero`][poly.standard.evaluation.polyvalzero]               | [`hermvalzero`][poly.hermite.hermvalzero]                 |
| Composition                | [`polycom`][poly.standard.evaluation.polycom]                       |                                                           |
|                            | [`polycom_naive`][poly.standard.evaluation.polycom_naive]           |                                                           |
|                            | [`polycom_iterative`][poly.standard.evaluation.polycom_iterative]   |                                                           |
|                            | [`polycom_horner`][poly.standard.evaluation.polycom_horner]         |                                                           |
| Shift                      | [`polyshift`][poly.standard.evaluation.polyshift]                   |                                                           |
| **Arithmetic**             |                                                                     |                                                           |
| Positive                   | [`polypos`][poly.standard.arithmetic.polypos]                       | [`hermpos`][poly.hermite.hermpos]                         |
| Negation                   | [`polyneg`][poly.standard.arithmetic.polyneg]                       | [`hermneg`][poly.hermite.hermneg]                         |
| Addition                   | [`polyadd`][poly.standard.arithmetic.polyadd]                       | [`hermadd`][poly.hermite.hermadd]                         |
| Basis addition             | [`polyaddc`][poly.standard.arithmetic.polyaddc]                     | [`hermaddc`][poly.hermite.hermaddc]                       |
| Subtraction                | [`polysub`][poly.standard.arithmetic.polysub]                       | [`hermsub`][poly.hermite.hermsub]                         |
| Scalar multiplication      | [`polyscalarmul`][poly.standard.arithmetic.polyscalarmul]           | [`hermscalarmul`][poly.hermite.hermscalarmul]             |
| Scalar true division       | [`polyscalartruediv`][poly.standard.arithmetic.polyscalartruediv]   | [`hermscalartruediv`][poly.hermite.hermscalartruediv]     |
| Scalar floor division      | [`polyscalarfloordiv`][poly.standard.arithmetic.polyscalarfloordiv] | [`hermscalarfloordiv`][poly.hermite.hermscalarfloordiv]   |
| Scalar mod                 | [`polyscalarmod`][poly.standard.arithmetic.polyscalarmod]           | [`hermscalarmod`][poly.hermite.hermscalarmod]             |
| Scalar divmod              | [`polyscalardivmod`][poly.standard.arithmetic.polyscalardivmod]     | [`hermscalardivmod`][poly.hermite.hermscalardivmod]       |
| Multiplication             | [`polymul`][poly.standard.arithmetic.polymul]                       | [`hermmul`][poly.hermite.hermmul]                         |
|                            | [`polymul_naive`][poly.standard.arithmetic.polymul_naive]           | [`hermmul_naive`][poly.hermite.hermmul_naive]             |
|                            | [`polymul_karatsuba`][poly.standard.arithmetic.polymul_karatsuba]   |                                                           |
| Multiplication by $x$      | [`polymulx`][poly.standard.arithmetic.polymulx]                     | [`hermmulx`][poly.hermite.hermmulx]                       |
| Exponentiation             | [`polypow`][poly.standard.arithmetic.polypow]                       | [`hermpow`][poly.hermite.hermpow]                         |
|                            | [`polypow_naive`][poly.standard.arithmetic.polypow_naive]           | [`hermpow_naive`][poly.hermite.hermpow_naive]             |
|                            | [`polypow_binary`][poly.standard.arithmetic.polypow_binary]         | [`hermpow_binary`][poly.hermite.hermpow_binary]           |
| Powers of polynomials      | [`polypows`][poly.standard.arithmetic.polypows]                     | [`hermpows`][poly.hermite.hermpows]                       |
| **Calculus**               |                                                                     |                                                           |
| Differentiation            | [`polyder`][poly.standard.calculus.polyder]                         | [`hermder`][poly.hermite.hermder]                         |
| Integration                | [`polyantider`][poly.standard.calculus.polyantider]                 | [`hermantider`][poly.hermite.hermantider]                 |
| **Conversion**             |                                                                     |                                                         |
| Sympification              | [`polysympify`][poly.standard.conversion.polysympify]               | [`hermsympify`][poly.hermite.hermsympify]                 |
| Unsympification            | [`polyunsympify`][poly.standard.conversion.polyunsympify]           |                                                           |

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
  - [ ] [`polymul_karatsuba`][poly.standard.polymul_karatsuba]
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
