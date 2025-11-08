# Standard basis

$$
    p(x)=\sum_{k=0}^na_kx^k \qquad P=(x_n)_{n\in\mathbb{N}_0}=(1, x, x^2, \dots)
$$

Prefixed by `poly...` (polynomial).

All functions accept **single exhaustible iterables**, if not stated otherwise.

They **return polynomials as `tuples`** of coefficients in **ascending order**.

Padding is done with `int(0)`.

If the result of a function is a polynomial, it is returned as a `tuple`.

Associative operations, like [`polyadd`][poly.functional.polyadd] and
[`polymul`][poly.functional.polymul], allow arbitrary many
arguments, even none.

The functions are **type-independent**. However, the data types used must
*support necessary scalar operations*. For instance, for polynomial addition,
components must be addable — this may include operations with padded integer
zeros. Empty operations return the zero polynomial (e.g. `polyadd()==polyzero`)
or try to infer the type from an argument (e.g.
`polyval(polyzero, x)==type(x)(0)`).

---

## Creation

::: poly.functional
    options:
      members:
        - polyzero
        - polyone
        - polyx
        - polymono
        - polymonos
        - polyrand
        - polyrandn
        - polyfromroots

## Utility

::: poly.functional
    options:
      members:
        - polyeq
        - polytrim
        - polyround
        - polydeg

## Evaluation

::: poly.functional
    options:
      members:
        - polyval

::: poly.functional
    options:
      heading_level: 4
      members:
        - polyval_naive
        - polyval_iterative
        - polyval_horner

::: poly.functional
    options:
      members:
        - polyvals
        - polyvalzero
        - polycom

::: poly.functional
    options:
      heading_level: 4
      members:
        - polycom_naive
        - polycom_iterative
        - polycom_horner

::: poly.functional
    options:
      members:
        - polycoms
        - polyshift

## Arithmetic

::: poly.functional
    options:
      members:
        - polypos
        - polyneg
        - polyadd
        - polyaddc
        - polysub
        - polyscalarmul
        - polyscalartruediv
        - polyscalarfloordiv
        - polyscalarmod
        - polyscalardivmod
        - polymul

::: poly.functional
    options:
      heading_level: 4
      members:
        - polymul_naive
        - polymul_karatsuba

::: poly.functional
    options:
      members:
        - polymulx
        - polypow

::: poly.functional
    options:
      heading_level: 4
      members:
        - polypow_naive
        - polypow_binary

## Calculus

::: poly.functional
    options:
      members:
        - polyder
        - polyantider

## Conversion

::: poly.functional
    options:
      members:
        - polysympify
        - polyunsympify
