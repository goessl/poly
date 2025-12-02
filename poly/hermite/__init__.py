r"""Polynomials in Hermite polynomial basis basis.

$$
    h(x) = \sum_{k=0}^nh_kH_k(x) \qquad H = (H_n)_{n\in\mathbb{N}_0} = (H_0, H_1, H_2, \dots)
$$

Polynomial results are returned in Hermite polynomial basis $H$, if not stated otherwise.

Integers are used where possible, for multiplications and for divisions like in [`hermantider`][poly.hermite.calculus.hermantider]. `Fraction`s are used where this was not possible, for example when a rational value has to be defined in [`hermx`][poly.hermite.hermx].
"""



from .creation import *
from .utility import *
from .evaluation import *
from .hilbert_space import *
from .arithmetic import *
from .calculus import *
from .conversion import *
