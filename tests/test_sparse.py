from poly import *
import numpy as np
from random import randint, random
from functools import reduce
from itertools import count, islice
from fractions import Fraction
from vector import vecrand
import sympy as sp
from sympy.abc import x as spx
import pytest



#creation
def test_polyszero():
    assert polyszero == {}

def test_polysone():
    assert polysone == {0:1}

def test_polysx():
    assert polysx == {1:1}

def test_polysmono():
    for n in range(10):
        for c in range(-10, +10):
            assert polysmono(n, c) == {n:c}
    assert polysmono(-1, 25) == polyszero

def test_polysmonos():
    for n, xn in enumerate(islice(polysmonos(), 10)):
        assert polysmono(n) == xn

def test_polysrand():
    for n in range(10):
        p = polysrand(n)
        assert isinstance(p, dict) and polysdeg(p)==n
        assert all(isinstance(ai, float) and 0<=ai<1 for ai in p.values())

def test_polysrandn():
    for _ in range(10):
        n = randint(0, 10)
        p = polysrandn(n)
        assert isinstance(p, dict) and polysdeg(p)==n
        assert all(isinstance(ai, float) for ai in p.values())
    assert polysrandn(-1, normed=False) == polyszero

def test_polysfromroots():
    for _ in range(100):
        n = randint(0, 10)
        r = vecrand(n)
        assert np.allclose(polystod(polysfromroots(*r)),
                                    polyfromroots(*r))
    assert polysfromroots() == polysone


#utility
def test_polysdeg():
    assert polysdeg(polyszero) == -1
    assert polysdeg(polysone) == 0
    assert polysdeg(polysx) == 1
    assert polysdeg(polysmono(3)) == 3

def test_polyseq():
    assert polyseq({}, {})
    assert polyseq({0:0}, {})
    assert polyseq({0:0, 1:0}, {0:0})
    assert polyseq({0:1}, {0:1, 1:0})
    assert not polyseq({0:1}, {})
    assert not polyseq({0:1, 1:2, 2:3}, {0:1, 1:2, 2:4})

def test_polystrim():
    assert polystrim({}) == {}
    assert polystrim({5:0}) == {}
    assert polystrim({2:0, 5:1}) == {5:1}

def test_polysround():
    assert polysround({}) == {}
    assert polysround({0:1.1, 1:2.2}) == {0:1, 1:2}
    assert polysround({0:1.12, 1:2.23}, ndigits=1) == {0:1.1, 1:2.2}


#evaluation
def test_polysval():
    #compare with numpy
    for _ in range(1000):
        p, x = polysrand(randint(1, 10)), random()
        assert np.isclose(polysval(p, x),
                          polyval(polystod(p), x))

def test_polysvalzero():
    assert polysvalzero(polyszero) == 0
    assert polysvalzero(polysone) == 1
    assert polysvalzero(polysx) == 0
    assert polysvalzero({0:5, 1:4, 2:3}) == 5

def test_polyscom():
    for _ in range(100):
        p, q = polysrand(randint(-1, 10)), polysrand(randint(-1, 10))
        assert np.allclose(polystod(polyscom(p, q)),
                           polycom(polystod(p), polystod(q)))

def test_polysshift():
    assert polysshift(polyszero, 5) == polyszero
    assert polysshift(polysone, 5) == polysone
    assert polysshift({0:1, 1:2, 2:3}, 5) == polyscom({0:1, 1:2, 2:3}, polysaddc(polysx, -5))

def test_polysscale():
    assert polysscale(polyszero, 5) == polyszero
    assert polysscale(polysone, 5) == polysone
    assert polysscale({0:1, 1:2, 2:3}, 5) == polyscom({0:1, 1:2, 2:3}, polysscalarmul(5, polysx))


#arithmetic
def test_polyspos():
    assert polyspos({0:+1, 1:-2, 2:+3}) == {0:+1, 1:-2, 2:+3}

def test_polysneg():
    assert polysneg({0:+1, 1:-2, 2:+3}) == {0:-1, 1:+2, 2:-3}

def test_polysadd():
    assert polysadd() == {}
    assert polysadd({0:1, 1:2}) == {0:1, 1:2}
    assert polysadd({0:1, 1:2, 2:3}, {0:4, 1:5}) == {0:5, 1:7, 2:3}

def test_polysaddc():
    assert polysaddc({}, 2, 3) == {3:2}
    assert polysaddc({0:1, 1:2}, 4, 5) == {0:1, 1:2, 5:4}
    assert polysaddc({0:1, 1:2, 2:3, 3:4, 4:5}, 5, 2) == {0:1, 1:2, 2:8, 3:4, 4:5}

def test_polyssub():
    assert polyssub({}, {}) == {}
    assert polyssub({0:1}, {}) == {0:1}
    assert polyssub({}, {0:1}) == {0:-1}
    assert polyssub({0:1, 1:2, 2:3}, {0:4, 1:5}) == {0:-3, 1:-3, 2:3}

def test_polysscalarmul():
    assert polysscalarmul(5, polyszero) == polyszero
    assert polysscalarmul(5, {0:1, 1:2, 2:3}) == {0:5, 1:10, 2:15}

def test_polysscalartruediv():
    assert polysscalartruediv({}, 1) == {}
    assert polysscalartruediv({5:4}, 2) == {5:2}

def test_polysscalarfloordiv():
    assert polysscalarfloordiv({}, 1) == {}
    assert polysscalarfloordiv({5:3}, 2) == {5:1}

def test_polysscalarmod():
    assert polysscalarmod({}, 2) == {}
    assert polysscalarmod({5:3}, 2) == {5:1}

def test_polysdivmod():
    p, a = {0:1, 1:2, 2:3}, 2
    assert polysscalardivmod(p, a) == (polysscalarfloordiv(p, a), polysscalarmod(p, a))

def test_polymul():
    for _ in range(1000):
        ps = [polysrand(randint(-1, 10)) for _ in range(randint(1, 4))]
        prediction = polysmul(*ps)
        actual = polymul(*map(polystod, ps))
        assert np.allclose(polystod(prediction), actual)

def test_polysmulx():
    for _ in range(1000):
        p = polysrand(randint(-1, 10))
        n = randint(0, 10)
        assert polyeq(polystod(polysmulx(p, n)), polymulx(polystod(p), n))

def test_polyspow():
    for method in {'naive', 'binary'}:
        for _ in range(1000):
            p, n = polysrand(randint(-1, 4)), randint(0, 3)
            assert np.allclose(polystod(polyspow(p, n, method)),
                                        polypow(polystod(p), n))

#calculus
def test_polysder():
    for _ in range(1000):
        p = polysrand(randint(-1, 10))
        k = randint(1, 10)
        assert np.allclose(polystod(polysder(p, k)),
                                    polyder(polystod(p), k))

def test_polysantider():
    for _ in range(1000):
        p = polysrand(randint(-1, 10))
        b, c = random(), random()
        assert np.allclose(polystod(polysantider(p, b=b, c=c)),
                                    polyantider(polystod(p), b=b, c=c))

#sympy
def test_polyssympify():
    assert polyssympify({0:1, 1:2, 2:3}) == sp.Poly(1+2*spx+3*spx**2)
    assert polyssympify(polyszero) == 0

def test_polysunsympify():
    assert polysunsympify(sp.Poly(1+2*spx+3*spx**2)) == {0:1, 1:2, 2:3}
