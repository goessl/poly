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
      heading_level: 3
      show_root_toc_entry: false
      show_symbol_type_heading: true
      members:
        - polyzero
        - polyone
        - polyx
        - polymono
        - polyrand
        - polyrandn
        - polyfromroots

---

## Utility

::: poly.functional
    options:
      heading_level: 3
      show_root_toc_entry: false
      show_symbol_type_heading: true
      members:
        - polyeq
        - polytrim
        - polyround
        - polydeg

---

## Evaluation

::: poly.functional
    options:
      heading_level: 3
      show_root_toc_entry: false
      show_symbol_type_heading: true
      members:
        - polyval
        - polyval_naive
        - polyval_iterative
        - polyval_horner
        - polyvalgen
        - polyvalzero
        - polycom
        - polycom_naive
        - polycom_iterative
        - polycom_horner
        - polycomgen
        - polyshift

---

## Arithmetic

::: poly.functional
    options:
      heading_level: 3
      show_root_toc_entry: false
      show_symbol_type_heading: true
      members:
        - polypos
        - polyneg
        - polyadd
        - polyaddc
        - polysub
        - polyscalarmul
        - polyscalartruediv
        - polyscalarfloordiv
        - polyfloordiv
        - polyscalarmod
        - polyscalardivmod
        - polymul
        - polymul_naive
        - polymul_karatsuba
        - polymulx
        - polypow
        - polypow_naive
        - polypow_binary

---

## Calculus

::: poly.functional
    options:
      heading_level: 3
      show_root_toc_entry: false
      show_symbol_type_heading: true
      members:
        - polyder
        - polyantider

---

## Conversion

::: poly.functional
    options:
      heading_level: 3
      show_root_toc_entry: false
      show_symbol_type_heading: true
      members:
        - polysympify
        - polyunsympify
