"""Coordinate two-jets and curvature for the numerical checks."""
from __future__ import annotations
import numpy as np

D = 6


class Jet:
    def __init__(self, value, first=None, second=None):
        self.v = np.asarray(value, dtype=float)
        if self.v.ndim == 0:
            self.v = self.v.reshape(1, 1)
        if self.v.ndim != 2:
            raise ValueError("all jet values are matrices, including column vectors")
        self.d = np.zeros((D,)+self.v.shape) if first is None else first
        self.dd = np.zeros((D, D)+self.v.shape) if second is None else second

    def __add__(self, other):
        other = asjet(other)
        return Jet(self.v+other.v, self.d+other.d, self.dd+other.dd)

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.d, -self.dd)

    def __sub__(self, other):
        return self + -asjet(other)

    def __rsub__(self, other):
        return asjet(other) + -self

    def __mul__(self, other):
        other = asjet(other)
        return Jet(
            self.v*other.v,
            self.d*other.v+self.v*other.d,
            self.dd*other.v+self.v*other.dd
            +self.d[:, None]*other.d[None, :]
            +self.d[None, :]*other.d[:, None],
        )

    __rmul__ = __mul__

    def __matmul__(self, other):
        other = asjet(other)
        cross = np.einsum("aik,bkj->abij", self.d, other.d)
        return Jet(
            self.v @ other.v,
            np.matmul(self.d, other.v)+np.matmul(self.v, other.d),
            np.matmul(self.dd, other.v)+np.matmul(self.v, other.dd)
            +cross+cross.swapaxes(0, 1),
        )

    def __pow__(self, power):
        if not isinstance(power, int) or power < 0:
            raise ValueError("only nonnegative integer scalar powers are used")
        answer = Jet(1.)
        for _ in range(power):
            answer = answer*self
        return answer

    @property
    def T(self):
        return Jet(self.v.T, self.d.transpose(0, 2, 1),
                   self.dd.transpose(0, 1, 3, 2))

    def entry(self, row, col=0):
        return Jet(self.v[row:row+1, col:col+1],
                   self.d[:, row:row+1, col:col+1],
                   self.dd[:, :, row:row+1, col:col+1])

    def inverse(self):
        P = np.linalg.inv(self.v)
        left = np.matmul(P, self.d)
        right = np.matmul(left, P)
        cross = np.matmul(left[:, None], right[None, :])
        return Jet(P, -right,
                   cross+cross.swapaxes(0, 1)-np.matmul(np.matmul(P, self.dd), P))


def asjet(value):
    return value if isinstance(value, Jet) else Jet(value)


def hstack(parts):
    parts = [asjet(p) for p in parts]
    return Jet(np.concatenate([p.v for p in parts], axis=1),
               np.concatenate([p.d for p in parts], axis=2),
               np.concatenate([p.dd for p in parts], axis=3))


def vstack(parts):
    parts = [asjet(p) for p in parts]
    return Jet(np.concatenate([p.v for p in parts], axis=0),
               np.concatenate([p.d for p in parts], axis=1),
               np.concatenate([p.dd for p in parts], axis=2))


def block(A, B, C, E):
    return vstack((hstack((A, B)), hstack((C, E))))


def sphere(p0, tangent, offset):
    first = np.zeros((D, 4, 1))
    second = np.zeros((D, D, 4, 1))
    jd = np.zeros((D, 4, 3))
    jdd = np.zeros((D, D, 4, 3))
    for a in range(3):
        first[offset+a, :, 0] = tangent[:, a]
        second[offset+a, offset+a, :, 0] = -p0
        jd[offset+a, :, a] = -p0
        for b in range(3):
            for c in range(3):
                jdd[offset+b, offset+c, :, a] = (
                    -tangent[:, a]*(b == c)
                    -tangent[:, b]*(a == c)
                    -tangent[:, c]*(a == b)
                )
    return Jet(p0[:, None], first, second), Jet(tangent, jd, jdd)


RIGHT_I = np.array([[0,-1,0,0],[1,0,0,0],[0,0,0,1],[0,0,-1,0]], dtype=float)
RIGHT_J = np.array([[0,0,-1,0],[0,0,0,-1],[1,0,0,0],[0,1,0,0]], dtype=float)
RIGHT_K = np.array([[0,0,0,-1],[0,0,1,0],[0,-1,0,0],[1,0,0,0]], dtype=float)


def left_frame(p):
    return hstack([Jet(R) @ p for R in (RIGHT_I, RIGHT_J, RIGHT_K)])


def curvature(g,dg,ddg):
    n=6; gi=np.linalg.inv(g)
    Gl=.5*(dg+dg.transpose(1,0,2)-dg.transpose(1,2,0))
    Gu=np.einsum('ijk,kl->ijl',Gl,gi)
    dGl=.5*(ddg+ddg.transpose(0,2,1,3)-ddg.transpose(0,2,3,1))
    dgi=-np.einsum('ij,ajk,kl->ail',gi,dg,gi)
    dGu=np.einsum('aijk,kl->aijl',dGl,gi)+np.einsum('ijk,akl->aijl',Gl,dgi)
    Ru=np.zeros((n,n,n,n))
    for i in range(n):
        for j in range(n):
            Ru[i,j]=dGu[i,j]-dGu[j,i]+np.einsum('km,ml->kl',Gu[j],Gu[i])-np.einsum('km,ml->kl',Gu[i],Gu[j])
    return np.einsum('ijkl,lm->ijkm',Ru,g)
