"""Exact algebraic checks for the proposed S^3 x S^3 construction.

All curvature calculations below are over QQ(z), with cos(theta)=(1-z^2)/(1+z^2)
and sin(theta)=2z/(1+z^2).  No floating-point arithmetic or input data files are
used.  The accompanying geometric argument is required to turn these algebraic
identities into the global positivity statement.
"""
from exact_engine import *
from exact_sparse import Linear, cubic_sparse
from exact_cubic import kbbjet, symprod
from exact_cross_checks import bo_Lgrad, L_on_flat
from exact_h0 import h0jet, at0, deriv0
from exact_diagonal import coef, arr0, h_k, base_hessian, extract_multilinear
import argparse


def assert_zero(value, name):
    a=np.asarray(value,dtype=object)
    if not all(x==0 for x in a.ravel()):
        raise AssertionError(name)


def first_family():
    start=time.time()
    c=(1-z*z)/(1+z*z);s=2*z/(1+z*z)
    r=ar([c,s,zero,zero]);p=eye(4)[0];q=r
    base=ExactBase(*g_jets(r));base.setup_plane()
    assert_zero(base.ddg-base.ddg.swapaxes(0,1)-eins('ijk,klm->ijlm',C,base.dg), 'frame commutator')
    assert_zero(base.R+base.R.swapaxes(0,1),'curvature skew first pair')
    assert_zero(base.R+base.R.swapaxes(2,3),'curvature skew second pair')
    assert_zero(base.R-base.R.transpose(2,3,0,1),'curvature pair symmetry')
    assert_zero(base.R+base.R.transpose(1,2,0,3)+base.R.transpose(2,0,1,3),'first Bianchi identity')
    assert base.D==Q(15,8)
    assert base.R[0,3,3,0]==0
    hj0=h0jet(r);linear0=Linear(base,hj0);data0=linear0.data()
    assert linear0.component(0,3,3,0)==0
    assert_zero(data0[0],'H0 connection on first-family torus')
    assert_zero(data0[1],'H0 curvature gradient on first family')
    print('PASS: background identities and first-family H0 cancellation',flush=True)

    # O-O plus the explicit K_OO correction.
    cs=coord_jets(p,q)
    height_lines=[Linear(base,const_height(cs,i,AP if i<4 else AQ)) for i in range(8)]
    data=[li.data() for li in height_lines]
    for li in height_lines:assert li.component(0,3,3,0)==0
    for ii in [0,1,4,5]:
        assert_zero(data[ii][1],'tangential-height first normal curvature gradient')
        assert_zero(data[ii][0][:,[1,2,4,5]],'tangential-height normal TT connection')
    PP=ar([[base.quad(data[i],data[j])for j in range(4)]for i in range(4)])
    QQm=ar([[base.quad(data[i+4],data[j+4])for j in range(4)]for i in range(4)])
    PQ=ar([[2*base.quad(data[i],data[j+4])for j in range(4)]for i in range(4)])
    qp=qmul(q,eye(4)[1])
    PQ+=(Q(-1,5)*np.outer(eye(4)[1],qp)+Q(-787,5880)*np.outer(q,p))/base.D
    x=s*s;den=31752000*(25+3*x)
    P=12632468+1767741*x+5670*x*x
    V=3941932+686019*x
    W=3521524+636555*x
    A=P/(2*den);Db=V*c/den;Jb=-W*s/den;d=-Q(787,11025)*c
    expected=zeros((4,4));expected[2,2]=expected[3,3]=A
    assert_zero(PP-expected,'O_p self coefficient')
    assert_zero(QQm-expected,'O_q self coefficient')
    ex=zeros((4,4));ex[0,0]=ex[1,1]=d;ex[2,2]=ex[3,3]=d+Db;ex[2,3]=Jb;ex[3,2]=-Jb
    assert_zero(PQ-ex,'corrected O_p O_q coefficient')
    assert Q(12632468-4627951,2*31752000*28)>Q(1,250)
    print('PASS: exact uniform first-family quadratic positivity margin',flush=True)

    bj=killing(p,q,Bco);cj=killing(p,q,Cco);bl=Linear(base,bj);bd=bl.data()
    cl=Linear(base,cj)
    assert bl.component(0,3,3,0)==0 and cl.component(0,3,3,0)==0
    assert_zero(cl.data()[0],'C TT connection')
    assert_zero(cl.data()[1],'C first normal curvature gradient')
    assert base.quad(bd,bd)==Q(128,845)*(c*c-s*s)
    phi=kbbjet(base,p,q)
    assert base.quad(bd,bd)+L_on_flat(phi,base)==0
    # Check the sparse first-variation evaluator against the full tensor evaluator.
    fullG,fullR=base.linear(bj)
    assert_zero(fullG-bl.A,'independent connection evaluators')
    for i,j,k,l in [(1,3,3,0),(2,3,3,0),(4,3,3,0),(5,3,3,0),
                    (0,1,3,0),(0,2,3,0),(0,4,3,0),(0,5,3,0)]:
        assert bl.component(i,j,k,l)==fullR[i,j,k,l]
    print('PASS: B-B coefficient and conformal quadratic correction',flush=True)

    for j,ang in [(1,c*c-s*s),(2,2*c*s)]:
        dj=killing(p,q,{(0,j):Q(2)})
        val=2*base.quad(bd,Linear(base,dj).data())
        corr=L_on_flat(symprod(base,bj,dj),base)
        assert val==Q(128,325)*ang
        assert corr==-Q(16,65)*ang
        assert val+Q(8,5)*corr==0
    print('PASS: B-C transverse-jet compatibility',flush=True)

    # Exact checks at four rational phases; the identities are linear in the
    # common phase after using the SO(4) axis-variation formulas.
    phases=[eye(4)[0],eye(4)[1],ar([Q(3,5),Q(4,5),zero,zero]),ar([Q(4,5),Q(3,5),zero,zero])]
    for p in phases:
        q=qmul(p,r);pu=qmul(p,eye(4)[1]);cs=coord_jets(p,q);bd=Linear(base,killing(p,q,Bco)).data()
        aa=-Q(64,455)*s*c;bb=Q(32,6825)*(41-30*s*s);dd=Q(352,6825)*s;ee=Q(1312,6825)*c
        targets=[aa*p[0]-bb*p[1],bb*p[0]+aa*p[1],dd*p[0]+ee*p[1],-ee*p[0]+dd*p[1]]
        for (a,b),expected in zip([(0,2),(0,3),(1,2),(1,3)],targets):
            dotjet=tuple(-y for y in const_height(cs,b if a==0 else b+4,AP if a==0 else AQ))
            val=2*base.quad(Linear(base,dotjet).data(),bd)
            vv,grad=bo_Lgrad(p,pu,q,c,a,b,base.D)
            assert val==expected
            assert vv==0 and val+grad==0
    print('PASS: B-O transverse-jet compatibility, exact rational phases',flush=True)

    print('Cubic setup at %.2f seconds' % (time.time()-start),flush=True)
    p=eye(4)[0];q=r;bj=killing(p,q,Bco);cj=killing(p,q,Cco)
    hj=tuple(a+b for a,b in zip(bj,cj));bc=symprod(base,bj,cj)
    kj=tuple(a+Q(8,5)*b for a,b in zip(kbbjet(base,p,q),bc))
    print('Cubic jets ready at %.2f seconds' % (time.time()-start),flush=True)
    q2,q3,eta=cubic_sparse(base,hj,kj)
    assert q2==0
    assert q3==Q(512,375)+Q(9728,21125)*(c*c-s*s)
    expected=ar([-Q(96,65),Q(32,65)*c/s,-(64-24*s*s)/65,Q(24,65)*s*c,
                  (64-24*s*s)/65,Q(24,65)*s*c,Q(96,65),Q(32,65)*c/s])
    assert_zero(eta-expected,'first-order minimizing plane')
    print('PASS: cubic coefficient 512/375 + (9728/21125) cos(2 theta)',flush=True)
    print('First-family checks complete in %.2f seconds' % (time.time()-start),flush=True)


def second_family():
    start=time.time();r=ar([2*z/(1+z*z),zero,(1-z*z)/(1+z*z),zero]);base=ExactBase(*g_jets(r))
    X=eye(6)[0];Y=eye(6)[3];N=eye(6)[[1,2,4,5]];base.X=X;base.Y=Y;base.N=N
    G=base.g;R=base.R;D=(X@G@X)*(Y@G@Y)-(X@G@Y)**2
    H8,rg=base_hessian(base)
    H9=zeros((9,9));H9[:8,:8]=ar([[at0(x)for x in row]for row in H8])
    H9[:8,8]=H9[8,:8]=ar([deriv0(x)for x in rg]);H9[8,8]=deriv0(R[0,3,3,0],2)
    lin=Linear(base,h0jet(r));ss,rr=lin.data();r9=np.r_[ar([at0(x)for x in rr]),deriv0(lin.component(0,3,3,0))]
    sh=ar([[at0(x)for x in row]for row in ss]);g=ar([[at0(x)for x in row]for row in G]);fixed=sh[1]@g@sh[1]-sh[0]@g@sh[2]
    inds=[0,1,2,3,4,5,6,8];H=H9[np.ix_(inds,inds)];ri=r9[inds];solution=-inv(H)@ri
    full=zeros(9);full[inds]=solution
    assert_zero(H9@full+r9,'second-family normal stationarity')
    assert at0(lin.component(0,3,3,0))==0
    assert at0(D)==Q(15,14) and fixed==Q(1,6)
    q=(fixed-ri@inv(H)@ri/2)/at0(D)
    assert q==Q(282877,2278125)
    assert q>Q(1,10)
    print('PASS: full base-and-plane second-family coefficient 282877/2278125',flush=True)
    print('Second-family checks complete in %.2f seconds' % (time.time()-start),flush=True)


def diagonal():
    start=time.time();p=eye(4)[0];X=eye(6)[0];Y=eye(6)[3];N=eye(6)[[1,2,4,5]]
    H10=zeros((10,10));r10=zeros(10)
    for ax in (2,3):
        r=ar([(1-z*z)/(1+z*z),zero,zero,zero]);r[ax]=2*z/(1+z*z)
        base=ExactBase(*g_jets(r));base.X=X;base.Y=Y;base.N=N
        h,k=h_k(base,p,r);lin=Linear(base,h);ss,rr=lin.data();H8,rg=base_hessian(base)
        H10[:8,:8]=arr0(H8);idx=ax+6
        H10[:8,idx]=H10[idx,:8]=ar([coef(v,1)for v in rg]);H10[idx,idx]=2*coef(base.R[0,3,3,0],2)
        r10[:8]=arr0(rr);r10[idx]=coef(lin.component(0,3,3,0),1)
    inds=[0,1,2,3,4,5,8,9];hh=H10[np.ix_(inds,inds)];rh=r10[inds]
    eta=zeros(10);eta[inds]=-inv(hh)@rh
    assert_zero(H10@eta+r10,'diagonal full normal stationarity')
    assert_zero(eta-ar([-Q(192,65),zero,-Q(64,65),zero,Q(64,65),zero,zero,zero,zero,-Q(16,65)]),'diagonal optimizer')
    aa,bb=eta[8],eta[9];den=1+(aa*aa+bb*bb)*z*z
    r=ar([(1-(aa*aa+bb*bb)*z*z)/den,zero,2*aa*z/den,2*bb*z/den])
    base=ExactBase(*g_jets(r));base.X=X;base.Y=Y;base.N=N
    h,k=h_k(base,p,r);lh=Linear(base,h);lk=Linear(base,k)
    A=eta[:4]@N;B=eta[4:8]@N;vs=[[X,A],[Y,B],[Y,B],[X,A]]
    def c0(i,j,k,l):return base.R[i,j,k,l]
    def c2(i,j,k,l):return lk.component(i,j,k,l)+lh.A[i,k]@base.g@lh.A[j,l]-lh.A[i,l]@base.g@lh.A[j,k]
    g=arr0(base.g);hh0=arr0(h[0]);sh=arr0(ar([lh.A[0,0],lh.A[0,3],lh.A[3,3]]));sk=arr0(ar([lk.A[0,0],lk.A[0,3],lk.A[3,3]]))
    fixed3=2*sh[1]@g@sk[1]-sh[0]@g@sk[2]-sh[2]@g@sk[0]-sh[1]@hh0@sh[1]+sh[0]@hh0@sh[2]
    n0=extract_multilinear(c0,vs,0)
    n1=extract_multilinear(c0,vs,1)+coef(lh.component(0,3,3,0),0)
    assert n0==0 and n1==0
    n2=extract_multilinear(c0,vs,2)+extract_multilinear(lh.component,vs,1)+coef(c2(0,3,3,0),0)
    n3=extract_multilinear(c0,vs,3)+extract_multilinear(lh.component,vs,2)+extract_multilinear(c2,vs,1)+fixed3
    D=g[0,0]*g[3,3]-g[0,3]**2
    assert n2==0 and n3/D==Q(115712,63375)
    print('PASS: diagonal full normal calculation, quadratic zero and cubic 115712/63375',flush=True)
    print('Diagonal checks complete in %.2f seconds' % (time.time()-start),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--part',choices=['all','first','second','diagonal'],default='all')
    args=parser.parse_args();start=time.time()
    if args.part in ['all','first']:first_family()
    if args.part in ['all','second']:second_family()
    if args.part in ['all','diagonal']:diagonal()
    print('ALL REQUESTED EXACT CHECKS PASSED (%.2f seconds)' % (time.time()-start),flush=True)
