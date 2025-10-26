# Hermite Polynomials

$$
    h(x) = \sum_{k=0}^nh_kH_k(x) \qquad H = (H_n)_{n\in\mathbb{N}_0} = (H_0, H_1, H_2, \dots)
$$

Polynomial results are returned in Hermite polynomial basis $H$, if not stated otherwise.

Integers are used where possible, for multiplications and for divisions like in [`hermantider`][poly.hermite_functional.hermantider]. `Fraction`s are used where this was not possible, for example when a rational value has to be defined in [`hermx`][poly.hermite_functional.hermx].

---

## Creation

::: poly.hermite_functional
    options:
      members:
        - H0
        - H1
        - H2
        - herm

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - herm_recursive
        - herm_iterative
        - herm_explicit

::: poly.hermite_functional
    options:
      members:
        - herms
        - hermzero
        - hermone
        - hermx
        - hermmono

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - hermmono_recursive
        - hermmono_iterative
        - hermmono_explicit

::: poly.hermite_functional
    options:
      members:
        - hermmonos
        - hermrand
        - hermrandn
        - hermfromroots

## Utility

::: poly.hermite_functional
    options:
      members:
        - hermeq
        - hermtrim
        - hermround
        - hermdeg

## Conversion

::: poly.hermite_functional
    options:
      members:
        - herm2poly

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - herm2poly_naive
        - herm2poly_iterative
        - herm2poly_clenshaw

::: poly.hermite_functional
    options:
      members:
        - poly2herm

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - poly2herm_naive
        - poly2herm_iterative
        - poly2herm_horner

## Evaluation

::: poly.hermite_functional
    options:
      members:
        - hermval

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - hermval_naive
        - hermval_iterative
        - hermval_clenshaw

::: poly.hermite_functional
    options:
      members:
        - hermvals
        - hermvalzero
        - hermvalzeros

## Arithmetic

::: poly.hermite_functional
    options:
      members:
        - hermpos
        - hermneg
        - hermadd
        - hermaddc
        - hermsub
        - hermscalarmul
        - hermscalartruediv
        - hermscalarfloordiv
        - hermscalarmod
        - hermscalardivmod
        - hermmul

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - hermmul_naive

::: poly.hermite_functional
    options:
      members:
        - hermmulx
        - hermpow

::: poly.hermite_functional
    options:
      heading_level: 4
      members:
        - hermpow_naive
        - hermpow_binary
    
::: poly.hermite_functional
    options:
      members:
        - hermmulpow

## Calculus

::: poly.hermite_functional
    options:
      members:
        - hermder
        - hermantider

## Sympy

::: poly.hermite_functional
    options:
      members:
        - hermsympify
