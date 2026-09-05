from exact_engine import *

def jetmul(A,B):
 a,da,dda=A;b,db,ddb=B
 c=a@b;dc=eins('aij,jk->aik',da,b)+eins('ij,ajk->aik',a,db)
 ddc=eins('abij,jk->abik',dda,b)+eins('ij,abjk->abik',a,ddb)+eins('aij,bjk->abik',da,db)+eins('bij,ajk->abik',da,db)
 return c,dc,ddc
def inversejet(base):
 w=base.ginv;dw=-eins('ij,ajk,kl->ail',w,base.dg,w)
 ddw=eins('ij,ajk,kl,blm,mn->abin',w,base.dg,w,base.dg,w)+eins('ij,bjk,kl,alm,mn->abin',w,base.dg,w,base.dg,w)-eins('ij,abjk,kl->abil',w,base.ddg,w)
 return w,dw,ddw
def symprod(base,A,B):
 j=jetmul(jetmul(A,inversejet(base)),B)
 return tuple(x+x.swapaxes(-1,-2)for x in j)
def scalar_times(s,A):
 f,df,ddf=s;a,da,dda=A
 return a*f,da*f+eins('a,ij->aij',df,a),dda*f+eins('a,bij->abij',df,da)+eins('b,aij->abij',df,da)+eins('ab,ij->abij',ddf,a)
def scalar_prod(A,B):
 a,da,dda=A;b,db,ddb=B
 return a*b,da*b+a*db,dda*b+a*ddb+np.outer(da,db)+np.outer(db,da)
def r0jet(p,q):
 cs=coord_jets(p,q);out=[zero,zeros(6),zeros((6,6))]
 for i in range(4):
  j=scalar_prod(tuple(c[...,i]for c in cs),tuple(c[...,i+4]for c in cs))
  out=[x+y for x,y in zip(out,j)]
 return tuple(out)
def kbbjet(base,p,q):
 sq=scalar_prod(r0jet(p,q),r0jet(p,q));phi=tuple(-Q(192,845)*a for a in sq);phi=(phi[0]+Q(96,845),phi[1],phi[2])
 return scalar_times(phi,(base.g,base.dg,base.ddg))
def curvature_series(base,jets,nmax=3,verbose=False):
 start=time.time()
 G=[base.g]+[j[0]for j in jets];dG=[base.dg]+[j[1]for j in jets];ddG=[base.ddg]+[j[2]for j in jets]
 Gam=[base.Gamma];dGam=[base.dGamma];Rup=[base.Rup];R=[base.R]
 for n in range(1,nmax+1):
  bn=B_lower(G[n],dG[n])if n<len(G)else zeros((6,6,6))
  dbn=B_lower(dG[n],ddG[n])if n<len(G)else zeros((6,6,6,6))
  for i in range(1,min(n,len(G)-1)+1):
   bn-=eins('lm,ijm->ijl',G[i],Gam[n-i]);dbn-=eins('alm,ijm->aijl',dG[i],Gam[n-i])+eins('lm,aijm->aijl',G[i],dGam[n-i])
  gn=eins('kl,ijl->ijk',base.ginv,bn)
  dgn=eins('kl,aijl->aijk',base.ginv,dbn-eins('alm,ijm->aijl',base.dg,gn))
  Gam.append(gn);dGam.append(dgn)
  rn=dgn-dgn.swapaxes(0,1)-eins('ijm,mkl->ijkl',C,gn)
  for i in range(n+1):rn+=eins('jkm,iml->ijkl',Gam[i],Gam[n-i])-eins('ikm,jml->ijkl',Gam[i],Gam[n-i])
  Rup.append(rn);rc=zeros((6,6,6,6))
  for i in range(min(n,len(G)-1)+1):rc+=eins('lm,ijkm->ijkl',G[i],Rup[n-i])
  R.append(rc)
  if verbose:print('series order',n,'done',time.time()-start,flush=True)
 return R,Gam

def moving_cubic(base,jets):
 Rs,Gs=curvature_series(base,jets);X,Y,N=base.X,base.Y,base.N
 rr=np.r_[2*eins('ijkl,ai,j,k,l->a',Rs[1],N,Y,Y,X),2*eins('ijkl,i,aj,k,l->a',Rs[1],X,N,Y,X)]
 eta=-base.Hinv@rr;A=eta[:4]@N;B=eta[4:]@N;XX=[X,A];YY=[Y,B]
 nums=zeros(4)
 for i in range(2):
  for j in range(2):
   for k in range(2):
    for l in range(2):
     m=i+j+k+l
     for n in range(m,4):nums[n]+=eins('ijkl,i,j,k,l->',Rs[n-m],XX[i],YY[j],YY[k],XX[l]).item()
 # Here Q0,Q1,Q2 vanish exactly, so cubic denominator is just D.
 return nums/base.D,eta,Rs,Gs
