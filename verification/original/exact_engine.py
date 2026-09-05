import numpy as np, time
from sympy.polys.fields import field
from sympy.polys.domains import QQ
from sympy.polys.fields import FracElement
for _name in ('__mul__','__rmul__','__add__','__radd__','__sub__','__rsub__','__truediv__','__rtruediv__'):
 _old=getattr(FracElement,_name)
 def _guard(self,other,_old=_old):
  if isinstance(other,np.ndarray):return NotImplemented
  return _old(self,other)
 setattr(FracElement,_name,_guard)
F,z=field('z',QQ)
zero=F.zero; one=F.one
Q=lambda a,b=1: F(QQ(a,b))
def ar(a):return np.array(a,dtype=object)
def zeros(shape):return np.full(shape,zero,dtype=object)
def eye(n):
 a=zeros((n,n))
 for i in range(n):a[i,i]=one
 return a
def inv(a):
 n=len(a);A=np.concatenate((a.copy(),eye(n)),axis=1)
 for i in range(n):
  j=next(j for j in range(i,n)if A[j,i]!=0)
  if j!=i:A[[i,j]]=A[[j,i]]
  A[i]=A[i]/A[i,i]
  for j in range(n):
   if j!=i and A[j,i]!=0:A[j]=A[j]-A[j,i]*A[i]
 return A[:,n:]
def eins(s,*a):return np.einsum(s,*a,optimize=True)
def skew(v):return ar([[zero,-v[2],v[1]],[v[2],zero,-v[0]],[-v[1],v[0],zero]])
I3=eye(3);J=ar([skew(x)for x in I3]);C=zeros((6,6,6))
for off in (0,3):
 for i in range(3):
  for j in range(3):C[off+i,off+j,off:off+3]=2*J[i][:,j]
def B_lower(g,dg):
 return Q(1,2)*(dg+dg.swapaxes(-3,-2)-np.moveaxis(dg,-3,-1)
  +eins('ijl,...lk->...ijk',C,g)-eins('jkl,...li->...ijk',C,g)+eins('kil,...lj->...ijk',C,g))
def algebraic_R(G):return eins('jkm,iml->ijkl',G,G)-eins('ikm,jml->ijkl',G,G)-eins('ijm,mkl->ijkl',C,G)
def adj(r):
 w=r[0];v=r[1:];return (w*w-v@v)*I3+2*np.outer(v,v)+2*w*skew(v)
def g_jets(r,lam=Q(1,3),tau=Q(3,10)):
 R=adj(r);a=(3*lam+3)/2;c=(3*lam-3)/2
 M=np.block([[a*I3,c*I3],[c*I3,a*I3]])
 W=inv(M)+tau*np.block([[I3,R],[R.T,I3]])
 dR=zeros((6,3,3));ddR=zeros((6,6,3,3))
 for i in range(3):
  dR[i]=-2*J[i]@R;dR[i+3]=2*R@J[i]
  for l in range(3):
   ddR[l,i]=4*J[i]@J[l]@R;ddR[l+3,i+3]=4*R@J[l]@J[i]
   ddR[l,i+3]=-4*J[l]@R@J[i];ddR[l+3,i]=-4*J[i]@R@J[l]
 dW=zeros((6,6,6));ddW=zeros((6,6,6,6))
 dW[:,:3,3:]=tau*dR;dW[:,3:,:3]=tau*dR.swapaxes(-1,-2)
 ddW[:,:,:3,3:]=tau*ddR;ddW[:,:,3:,:3]=tau*ddR.swapaxes(-1,-2)
 g=inv(W);dg=-eins('ij,ajk,kl->ail',g,dW,g)
 ddg=(eins('ij,ajk,kl,blm,mn->abin',g,dW,g,dW,g)+eins('ij,bjk,kl,alm,mn->abin',g,dW,g,dW,g)-eins('ij,abjk,kl->abil',g,ddW,g))
 return g,dg,ddg
class ExactBase:
 def __init__(self,g,dg,ddg):
  self.g=g;self.dg=dg;self.ddg=ddg;self.ginv=inv(g)
  B=B_lower(g,dg);dB=B_lower(dg,ddg)
  self.Gamma=eins('kl,ijl->ijk',self.ginv,B)
  self.dGamma=eins('kl,aijl->aijk',self.ginv,dB-eins('alm,ijm->aijl',dg,self.Gamma))
  self.Rup=self.dGamma-self.dGamma.swapaxes(0,1)+algebraic_R(self.Gamma)
  self.R=eins('lm,ijkm->ijkl',g,self.Rup)
 def linear(self,jet):
  h,dh,ddh=jet;B1=B_lower(h,dh);dB1=B_lower(dh,ddh)
  G1=eins('kl,ijl->ijk',self.ginv,B1-eins('lm,ijm->ijl',h,self.Gamma))
  dG1=eins('kl,aijl->aijk',self.ginv,dB1-eins('alm,ijm->aijl',dh,self.Gamma)-eins('lm,aijm->aijl',h,self.dGamma)-eins('alm,ijm->aijl',self.dg,G1))
  R1up=dG1-dG1.swapaxes(0,1)+eins('jkm,iml->ijkl',G1,self.Gamma)+eins('jkm,iml->ijkl',self.Gamma,G1)-eins('ikm,jml->ijkl',G1,self.Gamma)-eins('ikm,jml->ijkl',self.Gamma,G1)-eins('ijm,mkl->ijkl',C,G1)
  R1=eins('lm,ijkm->ijkl',self.g,R1up)+eins('lm,ijkm->ijkl',h,self.Rup)
  return G1,R1
 def setup_plane(self,axis=0):
  self.X=eye(6)[axis];self.Y=eye(6)[3+axis];inds=[i for i in range(6)if i not in(axis,3+axis)];self.N=eye(6)[inds]
  X=self.X;Y=self.Y;N=self.N;R=self.R
  HA=2*eins('ijkl,ai,j,k,bl->ab',R,N,Y,Y,N)
  HB=2*eins('ijkl,i,aj,bk,l->ab',R,X,N,N,X)
  HAB=2*(eins('ijkl,ai,bj,k,l->ab',R,N,N,Y,X)+eins('ijkl,ai,j,bk,l->ab',R,N,Y,N,X))
  H=np.block([[HA,HAB],[HAB.T,HB]]);self.H=(H+H.T)/2
  self.Hinv=inv(self.H)
  self.D=(X@self.g@X)*(Y@self.g@Y)-(X@self.g@Y)**2
 def linear_data(self,jet):
  G,R=self.linear(jet);X=self.X;Y=self.Y;N=self.N
  rr=np.r_[2*eins('ijkl,ai,j,k,l->a',R,N,Y,Y,X),2*eins('ijkl,i,aj,k,l->a',R,X,N,Y,X)]
  ss=ar([eins('ijk,i,j->k',G,X,X),eins('ijk,i,j->k',G,X,Y),eins('ijk,i,j->k',G,Y,Y)])
  return ss,rr,G,R
 def quad(self,A,B):
  sa,ra=A[:2];sb,rb=B[:2]
  return (sa[1]@self.g@sb[1]-(sa[0]@self.g@sb[2]+sa[2]@self.g@sb[0])/2-(ra@self.Hinv@rb)/2)/self.D

def qmul(p,q):return np.r_[p[0]*q[0]-p[1:]@q[1:],p[0]*q[1:]+q[0]*p[1:]+skew(p[1:])@q[1:]]
Rquat=ar([np.stack([qmul(e,eye(4)[i+1])for e in eye(4)],axis=1)for i in range(3)])
Lquat=ar([np.stack([qmul(eye(4)[i+1],e)for e in eye(4)],axis=1)for i in range(3)])
def coord_jets(p,q):
 f=np.r_[p,q];df=zeros((6,8));ddf=zeros((6,6,8))
 for off,quat,coff in [(0,p,0),(3,q,4)]:
  for i in range(3):
   df[off+i,coff:coff+4]=Rquat[i]@quat
   for l in range(3):ddf[off+l,off+i,coff:coff+4]=Rquat[i]@Rquat[l]@quat
 return f,df,ddf

def killing(p,q,co):
 W=zeros((12,6));dW=zeros((6,12,6));ddW=zeros((6,6,12,6))
 W[:3,:3]=I3;W[3:6,:3]=adj(p);W[6:9,3:]=I3;W[9:12,3:]=adj(q)
 for off,rowoff in [(0,3),(3,9)]:
  R=W[rowoff:rowoff+3,off:off+3]
  for i in range(3):
   dW[off+i,rowoff:rowoff+3,off:off+3]=2*R@J[i]
   for l in range(3):ddW[off+l,off+i,rowoff:rowoff+3,off:off+3]=4*R@J[l]@J[i]
 S=zeros((12,12))
 for (i,j),c in co.items():S[i,j]+=c/2;S[j,i]+=c/2
 H=W.T@S@W
 dH=eins('aij,ik,kl->ajl',dW,S,W)+eins('ij,ik,akl->ajl',W,S,dW)
 ddH=eins('abij,ik,kl->abjl',ddW,S,W)+eins('ij,ik,abkl->abjl',W,S,ddW)+eins('aij,ik,bkl->abjl',dW,S,dW)+eins('bij,ik,akl->abjl',dW,S,dW)
 return H,dH,ddH
Bco={(0,1):Q(1,2),(0,7):Q(3,26),(1,3):Q(1,2),(1,6):Q(-3,26),(1,9):Q(-3,26),(3,7):Q(3,26),(6,7):Q(-1,2),(7,9):Q(-1,2)}
Cco={(0,0):one}
def const_height(jets,i,A):return tuple(eins('...,ij->...ij',j[...,i],A)for j in jets)
AP=np.block([[I3,Q(-2,7)*I3],[Q(-2,7)*I3,zeros((3,3))]])
AQ=np.block([[zeros((3,3)),Q(-2,7)*I3],[Q(-2,7)*I3,I3]])
