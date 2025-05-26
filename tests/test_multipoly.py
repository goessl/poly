from poly import *
import numpy as np



def test_multipolyadd():
    for _ in range(100):
        pshape = np.random.randint(1, 4, np.random.randint(1, 4))
        qshape = np.random.randint(1, 4, np.random.randint(1, 4))
        p, q = np.random.randint(-100, +100, pshape), np.random.randint(-100, +100, qshape)
        r = multipolyadd(p, q)
        
        assert multipolysympify(r) == multipolysympify(p)+multipolysympify(q)

def test_multipolysub():
    for _ in range(100):
        pshape = np.random.randint(1, 4, np.random.randint(1, 4))
        qshape = np.random.randint(1, 4, np.random.randint(1, 4))
        p, q = np.random.randint(-100, +100, pshape), np.random.randint(-100, +100, qshape)
        r = multipolysub(p, q)
        
        assert multipolysympify(r) == multipolysympify(p)-multipolysympify(q)

def test_multipolymul():
    for _ in range(100):
        pshape = np.random.randint(1, 4, np.random.randint(1, 4))
        qshape = np.random.randint(1, 4, np.random.randint(1, 4))
        p, q = np.random.randint(-100, +100, pshape), np.random.randint(-100, +100, qshape)
        r = multipolymul(p, q)
        
        assert multipolysympify(r) == multipolysympify(p)*multipolysympify(q)
