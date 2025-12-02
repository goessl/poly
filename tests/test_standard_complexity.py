from poly import *
from random import randint, random
from collections import Counter
from operationcounter import OperationCounter, count_ops



#utility
def test_polyeq():
    for _ in range(100):
        n, m = randint(-1, 20), randint(-1, 20)
        p, q = polyrand(n), polyrand(m)
        p, q = OperationCounter.wrapCollection(p), OperationCounter.wrapCollection(q)
        with count_ops() as counts:
            polyeq(p, q)
        assert counts <= Counter({'eq':min(n, m)+1, 'bool':abs(n-m)})

def test_polytrim():
    for _ in range(100):
        n = randint(-1, 20)
        p = OperationCounter.wrapCollection(polyrand(n))
        with count_ops() as counts:
            polytrim(p)
        assert counts == Counter({'abs':n+1, 'gt':n+1})


#evaluation
def test_polyval_iterative():
    for _ in range(100):
        n = randint(-1, 20)
        p, x = polyrand(n), random()
        p, x = OperationCounter.wrapCollection(p), OperationCounter.wrap(x)
        with count_ops() as counts:
            polyval_iterative(p, x)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        #additions
        if n >= 0:
            assert counts['add'] == n
        if n <= 0:
            assert counts['add'] == 0
        #multiplications
        assert counts['mul'] == (2*n-1 if n>0 else 0)

def test_polyval_horner():
    for _ in range(100):
        n = randint(-1, 20)
        p, x = polyrand(n), random()
        p, x = OperationCounter.wrapCollection(p), OperationCounter.wrap(x)
        with count_ops() as counts:
            polyval_horner(p, x)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        #additions
        if n >= 0:
            assert counts == Counter({'add':n, 'mul':n})
        if n <= 0:
            assert counts == {}

def test_polycom_naive():
    for _ in range(100):
        n, m = randint(-1, 10), randint(-1, 10)
        p, q = OperationCounter.wrapCollection(polyrand(n)), OperationCounter.wrapCollection(polyrand(m))
        with count_ops() as counts:
            polycom_naive(p, q)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        assert counts['add'] == (round(n**3*m**2/6+n**2*m/2-n*m**2/6-n*m/2+n) if n>=0 and m>=0 else 0)
        assert counts['mul'] == (round(n**3*m**2/6+n**3*m/6+n**2*m-n*m**2/6-n*m/6+n**2/2+n/2) if n>=0 and m>=0 else 0)

def test_polycom_iterative():
    for _ in range(1000):
        n, m = randint(-1, 10), randint(-1, 10)
        p, q = OperationCounter.wrapCollection(polyrand(n)), OperationCounter.wrapCollection(polyrand(m))
        with count_ops() as counts:
            polycom_iterative(p, q)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        assert counts['add'] == (round(n**2*m**2/2+n**2*m/2-n*m**2/2-n*m/2+n) if n>=0 and m>=0 else 0)
        assert counts['mul'] == (round(n**2*m**2/2+n**2*m-n*m**2/2+n*m+2*n-m-1) if n>=1 and m>=0 else 0)

def test_polycom_horner():
    for _ in range(100):
        n, m = randint(-1, 10), randint(-1, 10)
        p, q = OperationCounter.wrapCollection(polyrand(n)), OperationCounter.wrapCollection(polyrand(m))
        with count_ops() as counts:
            polycom_horner(p, q)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'pos', 'add', 'mul'}
        assert counts['pos'] == (0 if m>=0 or n<0 else n)
        assert counts['add'] == (round(n**2*m**2/2-n*m**2/2+n) if n>=0 and m>=0 else 0)
        assert counts['mul'] == (round(n**2*m**2/2+n**2*m/2-n*m**2/2+n*m/2+n) if n>=0 and m>=0 else 0)

#arithmetic
def test_polymul_naive():
    for _ in range(100):
        n, m = randint(-1, 20), randint(-1, 20)
        p, q = polyrand(n), polyrand(m)
        p, q = OperationCounter.wrapCollection(p), OperationCounter.wrapCollection(q)
        with count_ops() as counts:
            polymul_naive(p, q)
        assert counts == Counter({'iadd':(n*m if n>=1 and m>=1 else 0), 'mul':(n+1)*(m+1)})

def test_polypow_naive():
    one = OperationCounter(1)
    for _ in range(100):
        n, k = randint(-1, 10), randint(0, 10)
        p = OperationCounter.wrapCollection(polyrand(n))
        with count_ops() as counts:
            polypow_naive(p, k, one=one)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        assert counts['add'] == (n**2*k*(k-1)//2 if n>=0 else 0)
        assert counts['mul'] == ((n*k+2)*(n+1)*(k-1)//2 if k>0 else 0)

def test_polypow_binary():
    one = OperationCounter(1)
    for _ in range(1000):
        n, k = randint(2, 20), randint(2, 20)
        p = OperationCounter.wrapCollection(polyrand(n))
        with count_ops() as counts:
            polypow_binary(p, k, one=one)
        counts = OperationCounter.grouped(counts)
        assert counts.keys() <= {'add', 'mul'}
        
        L = k.bit_length()
        w = k.bit_count()
        C = L + w - 1 #https://oeis.org/A056791
        B = 2*(2**L-1) + k-2**((k&-k).bit_length()-1) + sum(2**j * (k>>(j+1)).bit_count() for j in range(L) if ((k>>j)&0x01)==0x01)
        A = (4**L-1)//3 + sum(2**(i+j) for j in range(L) for i in range(j) if ((k>>j)&0x01)==((k>>i)&0x01)==0x01)
        assert counts['add'] == A*n**2
        assert counts['mul'] == A*n**2 + B*n + C


#calculus
def test_polyder():
    for _ in range(100):
        n = randint(-1, 20)
        k = randint(1, 10)
        p = OperationCounter.wrapCollection(polyrand(n))
        with count_ops() as counts:
            polyder(p, k)
        if k <= n:
            assert counts == Counter({'rmul':n-k+1})
        if k > n:
            assert counts == {}
