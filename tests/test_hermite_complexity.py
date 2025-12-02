from poly import *
from random import randint, random
from collections import Counter
from operationcounter import OperationCounter, count_ops


#conversion
def test_herm2poly_naive():
    for n in range(-1, 20):
        h = OperationCounter.wrapCollection(hermrand(n))
        with count_ops() as counts:
            herm2poly_naive(h)
        assert counts == Counter({'add': n*(n+1)//2,
                                  'mul': (n+1)*(n+2)//2})

def test_herm2poly_iterative():
    for n in range(-1, 20):
        h = OperationCounter.wrapCollection(hermrand(n))
        with count_ops() as counts:
            herm2poly_iterative(h)
        assert counts == Counter({'add': n*(n+1)//2,
                                  'mul': (n+1)*(n+2)//2})

def test_herm2poly_clenshaw():
    for n in range(-1, 20):
        h = OperationCounter.wrapCollection(hermrand(n))
        with count_ops() as counts:
            herm2poly_clenshaw(h)
        counts = OperationCounter.grouped(counts)
        assert counts == Counter({'add': n+1,
                                  'sub': n*(n-1)//2 if n>=0 else 0,
                                  'mul': n**2 if n>=0 else 0})

#Hilbert space
def test_hermabsqi():
    for n in range(-1, 20):
        h = OperationCounter.wrapCollection(hermrand(n))
        with count_ops() as counts:
            hermabsqi(h)
        assert counts == Counter({'add':(n if n>=0 else 0),
                                  'mul':2*(n+1)})

def test_hermdoti():
    for n in range(-1, 20):
        for m in range(-1, 20):
            g = OperationCounter.wrapCollection(hermrand(n))
            h = OperationCounter.wrapCollection(hermrand(m))
            with count_ops() as counts:
                hermdoti(g, h)
            assert counts == Counter({'add':(min(n, m) if n>=0 and m>=0 else 0),
                                      'mul':2*(min(n, m)+1)})

#arithmetic
def test_hermmul_naive():
    for _ in range(100):
        n, m = randint(-1, 10), randint(-1, 10)
        g, h = OperationCounter.wrapCollection(hermrand(n)), OperationCounter.wrapCollection(hermrand(m))
        with count_ops() as counts:
            hermmul_naive(g, h)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        N, M = min(n, m), max(n, m)
        assert counts['add'] == ((N+1)*(N+2)*(3*M+3-N)//6-(n+m+1) if n>=0 and m>=0 else 0)
        assert counts['mul'] == ((N+1)*(N+2)*(3*M+3-N)//3 if n>=0 and m>=0 else 0)

def test_hermmulx():
    for n in range(20):
        h = OperationCounter.wrapCollection(hermrand(n))
        with count_ops() as counts:
            hermmulx(h)
        counts = OperationCounter.grouped(counts)
        assert counts == Counter({'add': (n-1 if n>0 else 0),
                                  'mul': (n-1 if n>0 else 0),
                                  'truediv': n+1})
