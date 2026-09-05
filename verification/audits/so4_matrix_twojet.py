"""Floating-point second jets from exact differentiation identities.

There are no finite-difference steps in this module. Floating-point arithmetic
still prevents its output from being a proof certificate.
"""
from __future__ import annotations
import numpy as np
from so4_fd_probe import B4, HSTAR, star_wedge

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
STAR_B = np.einsum("ab,aij->bij", HSTAR, B4)
STAR_WEDGE = [[star_wedge(ei, ej) for ei in np.eye(4)] for ej in np.eye(4)]


def left_frame(p):
    return hstack([Jet(R) @ p for R in (RIGHT_I, RIGHT_J, RIGHT_K)])


def star_r_e(r, e_index):
    return sum((r.entry(i)*Jet(STAR_WEDGE[e_index][i]) for i in range(4)),
               Jet(np.zeros((4, 4))))


def cross_matrix(k):
    x, y, z = (k.entry(i) for i in range(3))
    return vstack((hstack((0., -z, y)), hstack((z, 0., -x)), hstack((-y, x, 0.))))


def joint_jets(p0, Tp, q0, Tq, matrices, labels, parameters=(1., 1., 1., 1.)):
    w1, w2, right_time, left_time = parameters
    p, Jp = sphere(p0, Tp, 0)
    q, Jq = sphere(q0, Tq, 3)
    gp, gq = Jp.T @ Jp, Jq.T @ Jq
    cp, cq = gp.inverse() @ Jp.T, gq.inverse() @ Jq.T
    round_g = block(w1*gp, np.zeros((3, 3)), np.zeros((3, 3)), w2*gq)
    K = hstack([vstack((cp @ Jet(B) @ p, cq @ Jet(B) @ q)) for B in B4])
    J = hstack([vstack((cp @ Jet(B) @ p, -(cq @ Jet(B) @ q))) for B in B4])
    # Right quaternion multiplication is anti-self-dual in this B4 basis.
    action_weight = .5*right_time*(np.eye(6)-HSTAR)+.5*left_time*(np.eye(6)+HSTAR)
    G = (round_g.inverse()+2*(K @ Jet(action_weight) @ K.T)).inverse()
    Lp, Lq = left_frame(p), left_frame(q)
    to_left = block(Lp.T @ Jp, np.zeros((3, 3)),
                    np.zeros((3, 3)), Lq.T @ Jq)
    Kf, Jf = G @ K, G @ J
    actions = []
    for r in (p+q, q-p):
        columns = []
        for a in range(4):
            A = star_r_e(r, a)
            columns.append(vstack((cp @ A @ p, cq @ A @ q)))
        actions.append(G @ hstack(columns))
    c = p.T @ q
    factors = ((6/5)*(2-c), -2*(1-c)*(2+c),
               (6/5)*(2+c), -2*(1+c)*(2-c))
    columns = []
    zero3 = np.zeros((3, 3))
    for matrix in matrices:
        T = Jet(matrix)
        Mp = (p.T @ T @ p)*np.eye(3) + Lp.T @ T @ Lp
        Mq = (q.T @ T @ q)*np.eye(3) + Lq.T @ T @ Lq
        for _, family, typ, power in labels:
            if family == "tt":
                if typ == 0:
                    raw = block(Mp, zero3, zero3, zero3)
                elif typ == 1:
                    raw = block(zero3, zero3, zero3, Mq)
                else:
                    cross = .5*(Mp+Mq)
                    raw = block(zero3, cross, cross, zero3)
                column = to_left.T @ raw @ to_left
            elif family.startswith("skew_"):
                k = Lp.T @ T @ p if family == "skew_p" else Lq.T @ T @ q
                skew = cross_matrix(k)
                raw = block(zero3, skew, -skew, zero3)
                column = (c**power)*(to_left.T @ raw @ to_left)
            else:
                v = T @ (p if family == "p" else q)
                Mv = hstack([Jet(B) @ v for B in STAR_B])
                source = actions[0 if typ in (0, 3) else 1] @ Mv
                target = Kf if typ in (0, 2) else Jf
                product = source @ target.T
                column = factors[typ]*(c**power)*.5*(product+product.T)
            columns.append(column)
    return (G.v, G.d, G.dd, np.stack([h.v for h in columns]),
            np.stack([h.d for h in columns], axis=1),
            np.stack([h.dd for h in columns], axis=2))


def odd_boundary_jets(p0, Tp, q0, Tq, vectors, degree, paired=True, common_weight=0,
                      central=False):
    """Equal-time seed, with intrinsic-first-null minus boundary values."""
    p, Jp = sphere(p0, Tp, 0)
    q, Jq = sphere(q0, Tq, 3)
    gp, gq = Jp.T @ Jp, Jq.T @ Jq
    cp, cq = gp.inverse() @ Jp.T, gq.inverse() @ Jq.T
    round_g = block(gp, np.zeros((3, 3)), np.zeros((3, 3)), gq)
    K = hstack([vstack((cp @ Jet(B) @ p, cq @ Jet(B) @ q)) for B in B4])
    J = hstack([vstack((cp @ Jet(B) @ p, -(cq @ Jet(B) @ q))) for B in B4])
    G = (round_g.inverse()+2*(K @ K.T)).inverse()
    if central:
        Lp, Lq = left_frame(p), left_frame(q)
        to_left = block(Lp.T @ Jp, np.zeros((3, 3)),
                        np.zeros((3, 3)), Lq.T @ Jq)
    Kf, Jf = G @ K, G @ J
    actions = []
    for r in (p+q, q-p):
        actions.append(G @ hstack([
            vstack((cp @ star_r_e(r, a) @ p, cq @ star_r_e(r, a) @ q))
            for a in range(4)
        ]))
    c = p.T @ q
    factors = ((6/5)*(2-c), -2*(1-c)*(2+c),
               (6/5)*(2+c), -2*(1+c)*(2-c))
    columns = []
    for v in vectors:
        Mv = Jet(np.stack([B @ v for B in STAR_B], axis=1))
        base = []
        for typ in range(4):
            source = actions[0 if typ in (0, 3) else 1] @ Mv
            target = Kf if typ in (0, 2) else Jf
            product = source @ target.T
            base.append(factors[typ]*.5*(product+product.T))
        local = []
        if paired:
            local.extend((base[0]+.6*base[1], base[2]+.6*base[3]))
        local.extend((c**power)*base[typ] for typ in range(4)
                     for power in range(1, degree+1))
        if central:
            for frame in (Lp, Lq):
                skew = cross_matrix(frame.T @ Jet(v[:, None]))
                raw = block(np.zeros((3, 3)), skew, -skew, np.zeros((3, 3)))
                tensor = to_left.T @ raw @ to_left
                local.extend((c**power)*tensor for power in range(degree+1))
        if common_weight == 0:
            multipliers = (Jet(1.),)
        elif common_weight == 2:
            multipliers = (p.entry(0)**2-p.entry(1)**2, 2*p.entry(0)*p.entry(1))
        else:
            raise ValueError("this audit uses only antipodal-odd tensor sectors")
        for multiplier in multipliers:
            columns.extend(multiplier*entry for entry in local)
    return (G.v, G.d, G.dd, np.stack([h.v for h in columns]),
            np.stack([h.d for h in columns], axis=1),
            np.stack([h.dd for h in columns], axis=2))
