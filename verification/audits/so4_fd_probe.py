#!/usr/bin/env python3
"""Finite-difference external probe for SO(4) four-slot Schur blocks.

This is a discovery tool only.  It is deliberately limited to small
parameter subspaces; all numerical signs must be rechecked symbolically or
with interval arithmetic before being used mathematically.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
import numpy as np

np.set_printoptions(precision=10, suppress=False)


PAIRS = [(i,j) for i in range(4) for j in range(i+1,4)]
B4=[]
for i,j in PAIRS:
    A=np.zeros((4,4)); A[i,j]=1.; A[j,i]=-1.; B4.append(A)
B4=np.asarray(B4)

def levi(i,j,k,l):
    if len({i,j,k,l})<4:return 0
    a=[i,j,k,l]; inv=sum(a[r]>a[s] for r in range(4) for s in range(r+1,4))
    return -1 if inv%2 else 1

HSTAR=np.zeros((6,6))
for b,W in enumerate(B4):
    st=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            st[i,j]=.5*sum(levi(i,j,k,l)*W[k,l] for k in range(4) for l in range(4))
    for a,(i,j) in enumerate(PAIRS): HSTAR[a,b]=st[i,j]

def wedge_coeff(a,b):
    W=np.outer(a,b)-np.outer(b,a)
    return np.asarray([W[i,j] for i,j in PAIRS])

def skew_from_coeff(c): return np.einsum('b,bij->ij',c,B4)
def star_wedge(a,b): return skew_from_coeff(HSTAR@wedge_coeff(a,b))

def sphere_chart(x,p0,T):
    d=np.sqrt(1.+x@x); v=p0+T@x
    return v/d, T/d-np.outer(v,x)/d**3

def vc(w,J): return np.linalg.solve(J.T@J,J.T@w)
def so(a,b): return .5*(np.outer(a,b)+np.outer(b,a))

def circle(theta):
    p=np.array([math.cos(theta),math.sin(theta),0.,0.])
    t=np.array([-math.sin(theta),math.cos(theta),0.,0.])
    T=np.stack((t,np.array([0.,0.,1.,0.]),np.array([0.,0.,0.,1.])),axis=1)
    return p,T

def base_hmap(coords,p0,Tp,q0,Tq,tcheeg=1.,scale=.5):
    p,Jp=sphere_chart(coords[:3],p0,Tp); q,Jq=sphere_chart(coords[3:],q0,Tq)
    G0=np.zeros((6,6)); G0[:3,:3]=Jp.T@Jp; G0[3:,3:]=Jq.T@Jq
    K=[]; J=[]
    for B in B4:
        K.append(np.r_[vc(B@p,Jp),vc(B@q,Jq)])
        J.append(np.r_[vc(B@p,Jp),-vc(B@q,Jq)])
    K=np.stack(K,axis=1); J=np.stack(J,axis=1)
    G=scale*np.linalg.inv(np.linalg.inv(G0)+tcheeg*K@K.T)
    c=p@q; rp=p+q; rm=q-p
    Kf=G@K; Jf=G@J
    H=np.zeros((96,6,6)); I4=np.eye(4)
    for a in range(4):
        Ap=star_wedge(rp,I4[a]); Am=star_wedge(rm,I4[a])
        kp=np.r_[vc(Ap@p,Jp),vc(Ap@q,Jq)]
        km=np.r_[vc(Am@p,Jp),vc(Am@q,Jq)]
        kpf=G@kp; kmf=G@km
        for b in range(6):
            k=a*6+b
            H[k]= (6./5.)*(2.-c)*so(kpf,Kf[:,b])
            H[24+k]= -2.*(1.-c)*(2.+c)*so(kmf,Jf[:,b])
            H[48+k]= (6./5.)*(2.+c)*so(kmf,Kf[:,b])
            H[72+k]= -2.*(1.+c)*(2.-c)*so(kpf,Jf[:,b])
    return G,H

def jets_fd(p0,Tp,q0,Tq,selected,hc=2e-4,tcheeg=1.,scale=.5):
    """Coordinate jets of the background and selected h-map columns."""
    z=np.zeros(6); G0,H0=base_hmap(z,p0,Tp,q0,Tq,tcheeg,scale)
    H0=H0[selected]
    n=6; r=len(selected)
    dg=np.zeros((n,6,6)); ddg=np.zeros((n,n,6,6))
    dh=np.zeros((n,r,6,6)); ddh=np.zeros((n,n,r,6,6))
    def ev(x): return base_hmap(x,p0,Tp,q0,Tq,tcheeg,scale)
    for i in range(n):
        e=np.zeros(n);e[i]=hc
        gp,hp=ev(e); gm,hm=ev(-e)
        dg[i]=(gp-gm)/(2*hc); ddg[i,i]=(gp-2*G0+gm)/(hc*hc)
        dh[i]=(hp[selected]-hm[selected])/(2*hc)
        ddh[i,i]=(hp[selected]-2*H0+hm[selected])/(hc*hc)
    for i in range(n):
        for j in range(i):
            ei=np.zeros(n);ej=np.zeros(n);ei[i]=hc;ej[j]=hc
            gpp,hpp=ev(ei+ej); gpm,hpm=ev(ei-ej)
            gmp,hmp=ev(-ei+ej); gmm,hmm=ev(-ei-ej)
            val=(gpp-gpm-gmp+gmm)/(4*hc*hc)
            ddg[i,j]=ddg[j,i]=val
            hv=(hpp[selected]-hpm[selected]-hmp[selected]+hmm[selected])/(4*hc*hc)
            ddh[i,j]=ddh[j,i]=hv
    return G0,dg,ddg,H0,dh,ddh

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

def kval(jets,lam,z):
    G0,dg0,ddg0,H0,dh0,ddh0=jets
    G=G0+np.einsum('p,pij->ij',lam,H0)
    dg=dg0+np.einsum('p,apij->aij',lam,dh0)
    ddg=ddg0+np.einsum('p,abpij->abij',lam,ddh0)
    R=curvature(G,dg,ddg)
    X=np.zeros(6);Y=np.zeros(6);X[0]=1.;Y[3]=1.
    ni=[1,2,4,5]; X[ni]+=z[:4];Y[ni]+=z[4:]
    num=np.einsum('i,j,k,l,ijkl',X,Y,Y,X,R)
    den=(X@G@X)*(Y@G@Y)-(X@G@Y)**2
    return num/den

def qmatrix(theta,phi,selected,hc=2e-4,hs=2e-4,hz=2e-4):
    p0,Tp=circle(theta);q0,Tq=circle(phi)
    jets=jets_fd(p0,Tp,q0,Tq,selected,hc)
    r=len(selected); zero=np.zeros(r+8); f0=kval(jets,np.zeros(r),np.zeros(8))
    L=np.zeros((r,r)); C=np.zeros((r,8)); H=np.zeros((8,8))
    for i in range(r):
        ei=np.zeros(r);ei[i]=hs
        L[i,i]=(kval(jets,ei,np.zeros(8))-2*f0+kval(jets,-ei,np.zeros(8)))/(hs*hs)
        for j in range(i):
            ej=np.zeros(r);ej[j]=hs
            L[i,j]=L[j,i]=(kval(jets,ei+ej,np.zeros(8))-kval(jets,ei-ej,np.zeros(8))-kval(jets,-ei+ej,np.zeros(8))+kval(jets,-ei-ej,np.zeros(8)))/(4*hs*hs)
    for i in range(r):
        ei=np.zeros(r);ei[i]=hs
        for a in range(8):
            za=np.zeros(8);za[a]=hz
            C[i,a]=(kval(jets,ei,za)-kval(jets,ei,-za)-kval(jets,-ei,za)+kval(jets,-ei,-za))/(4*hs*hz)
    for a in range(8):
        za=np.zeros(8);za[a]=hz
        H[a,a]=(kval(jets,np.zeros(r),za)-2*f0+kval(jets,np.zeros(r),-za))/(hz*hz)
        for b in range(a):
            zb=np.zeros(8);zb[b]=hz
            H[a,b]=H[b,a]=(kval(jets,np.zeros(r),za+zb)-kval(jets,np.zeros(r),za-zb)-kval(jets,np.zeros(r),-za+zb)+kval(jets,np.zeros(r),-za-zb))/(4*hz*hz)
    ew,ev=np.linalg.eigh(.5*(H+H.T)); tol=max(1e-9,np.max(np.abs(ew))*1e-7)
    inv=np.where(ew>tol,1./ew,0.); Hi=(ev*inv)@ev.T
    Q=.5*(L-C@Hi@C.T);Q=.5*(Q+Q.T)
    return Q,ew

def reduced_basis(raw_indices):
    """Return slot indices/coefficients for first-null P/Q and R/S pairs."""
    # effective coefficients at c=0 are 12/5 and -4, so ratio 3/5.
    out=[]; names=[]
    for k in raw_indices:
        if k<24:
            out.append([(k,1.),(24+k,.6)]); names.append(f'PQ{k}')
        else:
            kk=k-24
            out.append([(48+kk,1.),(72+kk,.6)]); names.append(f'RS{kk}')
    return out,names

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--n',type=int,default=3);ap.add_argument('--mode',default='small')
    args=ap.parse_args()
    # Canonical normal rows/column used by the first positive endpoint block.
    if args.mode=='small':
        raw=[12,18,36,42]  # two P/Q and two R/S first-null channels
    elif args.mode=='twistor':
        raw=list(range(0,24))+list(range(24,48))
    else:
        raw=[12,36]
    combos,names=reduced_basis(raw)
    # Build a selected raw map list and a linear transformation from combo
    # coefficients to raw slot coefficients.
    selected=[]
    for c in combos:
        for k,_ in c:
            if k not in selected:selected.append(k)
    # Compute Q in raw selected coordinates then pull back to combos.
    T=np.zeros((len(selected),len(combos)))
    for j,c in enumerate(combos):
        for k,a in c:T[selected.index(k),j]=a
    Qsum=np.zeros((len(combos),len(combos))); hspec=[]
    for i in range(args.n):
        th=2*math.pi*(i+.137)/args.n
        for j in range(args.n):
            ph=2*math.pi*(j+.419)/args.n
            Q,ew=qmatrix(th,ph,selected)
            Qsum += T.T@Q@T; hspec.append(ew.tolist())
            print('sample',i,j,flush=True)
    Qavg=Qsum/(args.n*args.n); vals=np.linalg.eigvalsh(.5*(Qavg+Qavg.T))
    out={'proof_status':'external finite-difference discovery only','mode':args.mode,'raw_indices':raw,'combo_names':names,'eigenvalues':vals.tolist(),'Q':Qavg.tolist(),'hessian':hspec}
    Path(args.output).parent.mkdir(parents=True,exist_ok=True);Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print('eigenvalues',vals,'pos',sum(vals>1e-6),'neg',sum(vals<-1e-6),flush=True)

if __name__=='__main__':main()

