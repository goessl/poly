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

| Operation                  | [poly](standard.md)                                                 | [sparse](sparse.md)                                                 |
| -------------------------- | ------------------------------------------------------------------- | ------------------------------------------------------------------- |
| **Creation**               |                                                                     |                                                                     |
| $0$                        | [`polyzero`][poly.standard.creation.polyzero]                       | [`polyszero`][poly.sparse.polyszero]                                |
| $1$                        | [`polyone`][poly.standard.creation.polyone]                         | [`polysone`][poly.sparse.creation.polysone]                         |
| $x$                        | [`polyx`][poly.standard.creation.polyx]                             | [`polysx`][poly.sparse.creation.polysx]                             |
| $x^n$                      | [`polymono`][poly.standard.creation.polymono]                       | [`polymono`][poly.sparse.creation.polysmono]                        |
| $(x^n)_{n\in\mathbb{N}_0}$ | [`polymonos`][poly.standard.creation.polymonos]                     | [`polysmonos`][poly.sparse.creation.polysmonos]                     |
| Uniform random             | [`polyrand`][poly.standard.creation.polyrand]                       | [`polysrand`][poly.sparse.creation.polysrandn]                      |
| Normal random              | [`polyrandn`][poly.standard.creation.polyrandn]                     | [`polysrandn`][poly.sparse.creation.polysrandn]                     |
| From roots                 | [`polyfromroots`][poly.standard.creation.polyfromroots]             | [`polysfromroots`][poly.sparse.creation.polysfromroots]             |
| **Utility**                |                                                                     |                                                                     |
| Degree                     | [`polydeg`][poly.standard.utility.polydeg]                          | [`polysdeg`][poly.sparse.utility.polysdeg]                          |
| Comparison                 | [`polyeq`][poly.standard.utility.polyeq]                            | [`polyseq`][poly.sparse.utility.polyseq]                            |
| Trimming                   | [`polytrim`][poly.standard.utility.polytrim]                        | [`polystrim`][poly.sparse.utility.polystrim]                        |
| **Conversion**             |                                                                     |                                                                     |
|                            |                                                                     | [`polystod`][poly.sparse.conversion.polystod]                       |
|                            |                                                                     | [`polydtos`][poly.sparse.conversion.polydtos]                       |
|                            | [`polysympify`][poly.standard.conversion.polysympify]               | [`polyssympify`][poly.sparse.conversion.polyssympify]               |
|                            | [`polyunsympify`][poly.standard.conversion.polyunsympify]           | [`polysunsympify`][poly.sparse.conversion.polysunsympify]           |
| **Evaluation**             |                                                                     |                                                                     |
| Evaluation                 | [`polyval`][poly.standard.evaluation.polyval]                       | [`polysval`][poly.sparse.evaluation.polysval]                       |
|                            | [`polyval_naive`][poly.standard.evaluation.polyval_naive]           |                                                                     |
|                            | [`polyval_iterative`][poly.standard.evaluation.polyval_iterative]   |                                                                     |
| Evaluation of basis        | [`polyvals`][poly.standard.evaluation.polyvals]                     |                                                                     |
| Evaluation at $x=0$        | [`polyvalzero`][poly.standard.evaluation.polyvalzero]               | [`polysvalzero`][poly.sparse.evaluation.polysvalzero]               |
| Composition                | [`polycom`][poly.standard.evaluation.polycom]                       | [`polyscom`][poly.sparse.evaluation.polyscom]                       |
|                            | [`polycom_naive`][poly.standard.evaluation.polycom_naive]           |                                                                     |
|                            | [`polycom_iterative`][poly.standard.evaluation.polycom_iterative]   |                                                                     |
|                            | [`polycom_horner`][poly.standard.evaluation.polycom_horner]         |                                                                     |
| Shift                      | [`polyshift`][poly.standard.evaluation.polyshift]                   | [`polysshift`][poly.sparse.evaluation.polysshift]                   |
| Scale                      | [`polyscale`][poly.standard.evaluation.polyscale]                   | [`polysscale`][poly.sparse.evaluation.polysscale]                   |
| **Arithmetic**             |                                                                     |                                                                     |
| Positive                   | [`polypos`][poly.standard.arithmetic.polypos]                       | [`polyspos`][poly.sparse.arithmetic.polyspos]                       |
| Negation                   | [`polyneg`][poly.standard.arithmetic.polyneg]                       | [`polysneg`][poly.sparse.arithmetic.polysneg]                       |
| Addition                   | [`polyadd`][poly.standard.arithmetic.polyadd]                       | [`polysadd`][poly.sparse.arithmetic.polysadd]                       |
| Basis addition             | [`polyaddc`][poly.standard.arithmetic.polyaddc]                     | [`polysaddc`][poly.sparse.arithmetic.polysaddc]                     |
| Subtraction                | [`polysub`][poly.standard.arithmetic.polysub]                       | [`polyssub`][poly.sparse.arithmetic.polyssub]                       |
| Scalar multiplication      | [`polyscalarmul`][poly.standard.arithmetic.polyscalarmul]           | [`polysscalarmul`][poly.sparse.arithmetic.polysscalarmul]           |
| Scalar true division       | [`polyscalartruediv`][poly.standard.arithmetic.polyscalartruediv]   | [`polysscalartruediv`][poly.sparse.arithmetic.polysscalartruediv]   |
| Scalar floor division      | [`polyscalarfloordiv`][poly.standard.arithmetic.polyscalarfloordiv] | [`polysscalarfloordiv`][poly.sparse.arithmetic.polysscalarfloordiv] |
| Scalar mod                 | [`polyscalarmod`][poly.standard.arithmetic.polyscalarmod]           | [`polysscalarmod`][poly.sparse.arithmetic.polysscalarmod]           |
| Scalar divmod              | [`polyscalardivmod`][poly.standard.arithmetic.polyscalardivmod]     | [`polysscalardivmod`][poly.sparse.arithmetic.polysscalardivmod]     |
| Multiplication             | [`polymul`][poly.standard.arithmetic.polymul]                       | [`polysmul`][poly.sparse.arithmetic.polysmul]                       |
|                            | [`polymul_naive`][poly.standard.arithmetic.polymul_naive]           | [`polysmul_naive`][poly.sparse.arithmetic.polysmul_naive]           |
|                            | [`polymul_karatsuba`][poly.standard.arithmetic.polymul_karatsuba]   |                                                                     |
| Multiplication by $x$      | [`polymulx`][poly.standard.arithmetic.polymulx]                     | [`polysmulx`][poly.sparse.arithmetic.polysmulx]                     |
| Exponentiation             | [`polypow`][poly.standard.arithmetic.polypow]                       | [`polyspow`][poly.sparse.arithmetic.polyspow]                       |
|                            | [`polypow_naive`][poly.standard.arithmetic.polypow_naive]           | [`polyspow_naive`][poly.sparse.arithmetic.polyspow_naive]           |
|                            | [`polypow_binary`][poly.standard.arithmetic.polypow_binary]         | [`polyspow_binary`][poly.sparse.arithmetic.polyspow_binary]         |
| Powers of polynomials      | [`polypows`][poly.standard.arithmetic.polypows]                     | [`polyspows`][poly.sparse.arithmetic.polyspows]                     |
| **Calculus**               |                                                                     |                                                                     |
| Differentiation            | [`polyder`][poly.standard.calculus.polyder]                         | [`polysder`][poly.sparse.calculus.polysder]                         |
| Integration                | [`polyantider`][poly.standard.calculus.polyantider]                 | [`polysantider`][poly.sparse.calculus.polysantider]                 |

## Design

### Coefficient order

Coefficients are stored in ascending order like every sane person would do.

This ways the coefficient indices correspond to the monomial exponent (or the basis index incase of Hermite polynomials).

[The newer `numpy.polynomial` does it this way, the older `numpy.poly1d` does it in reverse](https://numpy.org/doc/stable/reference/routines.polynomials.html#transition-guide). [`sympy` also does it the wrong way](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.Poly.coeffs). 

## Roadmap

- Modules
  - [x] Hermite polynomials module
  - [x] Sparse module
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
