# Hermite Polynomials

$$
    h(x) = \sum_{k=0}^nh_kH_k(x) \qquad H = (H_n)_{n\in\mathbb{N}_0} = (H_0, H_1, H_2, \dots)
$$

Polynomial results are returned in Hermite polynomial basis $H$, if not stated otherwise.

Integers are used where possible, for multiplications and for divisions like in [`hermantider`][poly.hermite.hermantider]. `Fraction`s are used where this was not possible, for example when a rational value has to be defined in [`hermx`][poly.hermite.hermx].

---

## Creation

::: poly.hermite
    options:
      members:
        - H0
        - H1
        - H2
        - herm

::: poly.hermite
    options:
      heading_level: 4
      members:
        - herm_recursive
        - herm_iterative
        - herm_explicit
        - herms

::: poly.hermite
    options:
      members:
        - hermzero
        - hermone
        - hermx
        - hermmono

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermmono_recursive
        - hermmono_iterative
        - hermmono_explicit
        - hermmonos

::: poly.hermite
    options:
      members:
        - hermrand
        - hermrandn
        - hermfromroots

## Utility

::: poly.hermite
    options:
      members:
        - hermeq
        - hermtrim
        - hermround
        - hermdeg

## Hilbert space

::: poly.hermite
    options:
      members:
        - hermweight

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermweights
        - hermweighti
        - hermweightis

::: poly.hermite
    options:
      members:
        - hermabs

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermabsq
        - hermabsqi

::: poly.hermite
    options:
      members:
        - hermdot

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermdoti

## Conversion

::: poly.hermite
    options:
      members:
        - herm2poly

::: poly.hermite
    options:
      heading_level: 4
      members:
        - herm2poly_naive
        - herm2poly_iterative
        - herm2poly_clenshaw

::: poly.hermite
    options:
      members:
        - poly2herm

::: poly.hermite
    options:
      heading_level: 4
      members:
        - poly2herm_naive
        - poly2herm_iterative
        - poly2herm_horner

## Evaluation

::: poly.hermite
    options:
      members:
        - hermval

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermval_naive
        - hermval_iterative
        - hermval_clenshaw

::: poly.hermite
    options:
      members:
        - hermvals
        - hermvalzero
        - hermvalzeros

## Arithmetic

::: poly.hermite
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

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermmul_naive
        - hermmulx
        - hermmuln

::: poly.hermite
    options:
      members:
        - hermpow

::: poly.hermite
    options:
      heading_level: 4
      members:
        - hermpow_naive
        - hermpow_binary
        - hermmulpow

## Calculus

::: poly.hermite
    options:
      members:
        - hermder
        - hermantider

## Sympy

::: poly.hermite
    options:
      members:
        - hermsympify
