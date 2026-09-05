"""Supplementary exact audits independent of the coefficient-table interface."""
from exact_engine import *
from exact_cubic import scalar_times, scalar_prod, r0jet
from exact_cross_checks import boco,pairs6
from exact_sparse import Linear

def addj(*jj):return tuple(sum((j[k] for j in jj),zeros(jj[0][k].shape)) for k in range(3))
def scalej(a,j):return tuple(a*x for x in j)
def constant_tensor(co,p,q):return killing(p,q,co)

def K1_compact(p,q):
    cs=coord_jets(p,q);coords=[tuple(j[...,i]for j in cs)for i in range(8)]
    rr=r0jet(p,q)
    # A and B are expanded from the four compact coframe expressions in the proof.
    def loc(off, kind, xv):
        # Local ordering left i,j,k, right i,j,k.
        ix=[off,off+1,off+2,off+3,off+4,off+5]
        if kind=='A':
            co0={(ix[a],ix[b]):v for (a,b),v in { (0,1):one,(0,4):-one,(1,3):one,(3,4):-one }.items()}
            co1={(ix[a],ix[b]):v for (a,b),v in { (0,2):Q(13,7),(0,5):-one,(2,3):Q(13,7),(3,5):-one }.items()}
            i2=3
        else:
            co0={(ix[a],ix[b]):v for (a,b),v in { (0,2):Q(13,7),(0,5):one,(2,3):Q(13,7),(3,5):one }.items()}
            co1={(ix[a],ix[b]):-one for a,b in [(0,1),(0,4),(1,3),(3,4)]}
            i2=2
        co2={(ix[a],ix[b]):Q(-4,7) for a,b in [(0,0),(0,3),(3,3)]}
        return addj(scalar_times(xv[0],killing(p,q,co0)),scalar_times(xv[1],killing(p,q,co1)),scalar_times(xv[i2],killing(p,q,co2)))
    return addj(scalej(Q(11,325),addj(loc(6,'A',coords[:4]),loc(0,'B',coords[4:]))),
                scalej(Q(6,65),scalar_times(rr,addj(loc(0,'B',coords[:4]),loc(6,'A',coords[4:])))))

def K1_table(p,q):
    cs=coord_jets(p,q);sc=[tuple(j[...,i]for j in cs)for i in range(8)]
    rr=r0jet(p,q);sc+= [scalar_prod(rr,j)for j in sc]
    jj=[]
    for k,row in enumerate(boco):
        off=6 if k<4 or k>=12 else 0
        co={(off+i,off+j):a for a,(i,j)in zip(row,pairs6)if a}
        jj.append(scalar_times(sc[k],killing(p,q,co)))
    return addj(*jj)

def checkzero(a,label):
    assert all(x==0 for x in np.asarray(a,dtype=object).ravel()),label

def check_K1_formal():
    """Compare the compact formula and table as polynomials in 21 variables."""
    from sympy import Poly, Rational, symbols
    p = symbols('p0:4')
    q = symbols('q0:4')
    r0 = symbols('r0')
    coframes = symbols('a0:6') + symbols('b0:6')

    def local(forms, x, kind):
        l1, l2, l3, h1, h2, h3 = forms
        u = l1 + h1
        t = l1*l1 + l1*h1 + h1*h1
        if kind == 'A':
            return u*(x[0]*(l2-h2) + x[1]*(Rational(13,7)*l3-h3)) - Rational(4,7)*x[3]*t
        return u*(x[0]*(Rational(13,7)*l3+h3) - x[1]*(l2+h2)) - Rational(4,7)*x[2]*t

    compact = Rational(11,325)*(local(coframes[6:],p,'A') + local(coframes[:6],q,'B'))
    compact += Rational(6,65)*r0*(local(coframes[:6],p,'B') + local(coframes[6:],q,'A'))
    scalars = p + q + tuple(r0*x for x in p+q)
    table = 0
    for k, row in enumerate(boco):
        offset = 6 if k < 4 or k >= 12 else 0
        for coefficient, (i, j) in zip(row, pairs6):
            table += scalars[k]*Rational(str(coefficient))*coframes[offset+i]*coframes[offset+j]
    variables = p + q + (r0,) + coframes
    left = Poly(compact, *variables, domain='QQ')
    right = Poly(table, *variables, domain='QQ')
    if left != right:
        raise AssertionError('compact K1 and coefficient table differ as formal polynomials')
    return len(left.terms())

if __name__=='__main__':
    monomials = check_K1_formal()
    print(f'PASS: compact K1 equals its coefficient table as a formal polynomial ({monomials} monomials)',flush=True)
    p=ar([Q(1,2),Q(1,2),Q(1,2),Q(1,2)])
    q=ar([(1-z*z)/(1+z*z),2*z/(1+z*z),zero,zero])
    for a,b in zip(K1_compact(p,q),K1_table(p,q)):checkzero(a-b,'compact K1 versus coefficient table')
    print('PASS: compact K1 and table derivatives agree on the specified test curve',flush=True)
    gg=3*eye(6);base=ExactBase(gg,zeros((6,6,6)),zeros((6,6,6,6)))
    for i in range(6):
        for j in range(6):
            assert base.R[i,j,j,i] == (Q(3) if i!=j and i//3==j//3 else zero)
    print('PASS: curvature engine recovers product sectional curvature 1/3 and mixed zero',flush=True)
    # Independent full-connection recurrence at a rational angle.  This does
    # not use the Gram-term formula in the sparse cubic evaluator.
    import exact_cubic
    from exact_cubic import curvature_series,kbbjet,symprod
    from exact_sparse import cubic_sparse
    exact_cubic.t0=time.time()
    p=eye(4)[0];q=ar([Q(3,5),Q(4,5),zero,zero]);base=ExactBase(*g_jets(q));base.setup_plane()
    bj=killing(p,q,Bco);cj=killing(p,q,Cco)
    hj=addj(bj,cj);kj=addj(kbbjet(base,p,q),scalej(Q(8,5),symprod(base,bj,cj)))
    Rs,Gs=curvature_series(base,[hj,kj],nmax=3)
    q2,q3,eta=cubic_sparse(base,hj,kj)
    X,Y,N=base.X,base.Y,base.N;A=eta[:4]@N;B=eta[4:]@N
    nums=zeros(4)
    import itertools
    for ii,jj,kk,ll in itertools.product(range(2),repeat=4):
        deg=ii+jj+kk+ll
        for n in range(deg,4):
            nums[n]+=eins('ijkl,i,j,k,l->',Rs[n-deg],[X,A][ii],[Y,B][jj],[Y,B][kk],[X,A][ll]).item()
    assert nums[0]==0 and nums[1]==0 and nums[2]==0
    assert nums[3]/base.D==q3==Q(512,375)+Q(9728,21125)*Q(-7,25)
    print('PASS: independent full-connection curvature recurrence reproduces the cubic exactly',flush=True)
