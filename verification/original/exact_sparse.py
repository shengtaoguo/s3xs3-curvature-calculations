from exact_engine import *
from functools import lru_cache
class Linear:
 def __init__(self,base,jet):
  self.base=base;self.h,self.dh,self.ddh=jet
  b=B_lower(self.h,self.dh)-eins('lm,ijm->ijl',self.h,base.Gamma)
  self.A=eins('kl,ijl->ijk',base.ginv,b)
 @lru_cache(None)
 def derivative(self,i,j,k):
  base=self.base;h,dh,ddh=self.h,self.dh,self.ddh
  db=(ddh[i,j,k,:]+ddh[i,k,j,:]-ddh[i,:,j,k]
    +C[j,k]@dh[i]-np.einsum('lm,m->l',C[k],dh[i,:,j],optimize=True)+np.einsum('lm,m->l',C[:,j],dh[i,:,k],optimize=True))/2
  return base.ginv@(db-dh[i]@base.Gamma[j,k]-h@base.dGamma[i,j,k]-base.dg[i]@self.A[j,k])
 @lru_cache(None)
 def component(self,i,j,k,l):
  G=self.base.Gamma;A=self.A
  vec=self.derivative(i,j,k)-self.derivative(j,i,k)
  for m in range(6):vec+=A[j,k,m]*G[i,m]+G[j,k,m]*A[i,m]-A[i,k,m]*G[j,m]-G[i,k,m]*A[j,m]-C[i,j,m]*A[m,k]
  return self.base.g[l]@vec+self.h[l]@self.base.Rup[i,j,k]
 def data(self):
  axis=next(i for i,x in enumerate(self.base.X)if x);yi=axis+3
  inds=[i for i in range(6)if i not in(axis,yi)]
  rr=ar([*[2*self.component(i,yi,yi,axis)for i in inds],*[2*self.component(axis,i,yi,axis)for i in inds]])
  ss=ar([self.A[axis,axis],self.A[axis,yi],self.A[yi,yi]])
  return ss,rr

def nonzero(v):return[(i,x)for i,x in enumerate(v)if x]
def multilinear(component,X,Y,Z,W):
 out=zero
 for i,x in nonzero(X):
  for j,y in nonzero(Y):
   for k,zv in nonzero(Z):
    for l,w in nonzero(W):out+=component(i,j,k,l)*x*y*zv*w
 return out

def cubic_sparse(base,H,K):
 lh=Linear(base,H);lk=Linear(base,K);sh,rh=lh.data();eta=-base.Hinv@rh
 X,Y,N=base.X,base.Y,base.N;A=eta[:4]@N;B=eta[4:]@N
 sk,_=lk.data();h=H[0];g=base.g
 c3=(2*sh[1]@g@sk[1]-sh[0]@g@sk[2]-sh[2]@g@sk[0]-sh[1]@h@sh[1]+sh[0]@h@sh[2])
 def r2(i,j,k,l):
  aa=lh.A
  return lk.component(i,j,k,l)+aa[i,k]@g@aa[j,l]-aa[i,l]@g@aa[j,k]
 c3+=2*multilinear(r2,A,Y,Y,X)+2*multilinear(r2,X,B,Y,X)
 XX=[X,A];YY=[Y,B]
 for i in range(2):
  for j in range(2):
   for k in range(2):
    for l in range(2):
     deg=i+j+k+l
     if deg==2:c3+=multilinear(lh.component,XX[i],YY[j],YY[k],XX[l])
     if deg==3:c3+=eins('ijkl,i,j,k,l->',base.R,XX[i],YY[j],YY[k],XX[l]).item()
 c2=(sk[1]@g*0).sum() # allocate scalar
 c2=base.quad(lh.data(),lh.data())+multilinear(lk.component,X,Y,Y,X)/base.D
 return c2,c3/base.D,eta
