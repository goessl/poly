import numpy as np
from vector import vecnpadd

from ._hermfunctions import herm



__all__ = ['hermnp2poly', 'hermnp2polyM',
           'polynp2hermM', 'hermnpder', 'hermnpderM']



def hermnp2poly(h):
    return vecnpadd(*(hi.T*herm(i) for i, hi in enumerate(h.T[:,np.newaxis])))

def hermnp2polyM(deg):
    M = np.zeros((deg+1, deg+1), dtype=object)
    #for i in range(M.shape[1]):
    #    M[:i+1, i] = herm(i) #checks for overflow: at deg=25
    #return M
    M[0, 0] = 1
    #diagonal: powers of 2
    for i in range(1, deg+1):
        M[i, i] = 2 * M[i-1, i-1]
    #first row: Hermite numbers
    for n in range(2, deg+1, 2):
        M[0, n] = -2*(n-1) * M[0, n-2]
    #inner upper triangle: recursion relation
    for c in range(3, M.shape[1]):
        for r in range(int(bool(c%2)), c-1, 2):
            M[r, c] = 2*M[r-1, c-1] - (r+1)*M[r+1, c-1]
    return M

def polynp2hermM(deg):
    M = np.full((deg+1, deg+1), Fraction(0))
    M[0, 0] = Fraction(1)
    #diagonal: 2^(-n)
    for i in range(1, deg+1):
        M[i, i] = M[i-1, i-1] / 2
    #first row: n!/((n/2)!2^n)
    for n in range(2, deg+1, 2):
        #M[0, n] = Fraction(factorial(n), factorial(n//2)*2**n)
        M[0, n] = M[0, n-2] * (n-1) / 2
    #inner upper triangle: recursion relation
    for c in range(3, M.shape[1]):
        for r in range(int(bool(c%2)), c-1, 2):
            M[r, c] = M[r-1, c-1]/2 + (r+1)*M[r+1, c-1]
    return M




def hermnpder(h):
    h = np.asarray(h)
    return h[...,1:] * (np.arange(2, 2*h.shape[-1], 2))

def hermnpderM(deg):
    D = np.zeros((deg, deg+1), dtype=np.uint64)
    np.fill_diagonal(D[:,1:], np.arange(2, 2*D.shape[0]+1, 2))
    return D
