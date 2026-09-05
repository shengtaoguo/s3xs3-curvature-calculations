from exact_engine import *
from exact_sparse import Linear,nonzero
from exact_cubic import kbbjet,symprod

def coef(f,n):
 if not isinstance(f,FracElement):return F(f)if n==0 else zero
 a=[f.numer.get((i,),QQ.zero)for i in range(n+1)];b=[f.denom.get((i,),QQ.zero)for i in range(n+1)];out=[]
 for k in range(n+1):out.append((a[k]-sum(b[j]*out[k-j]for j in range(1,k+1)))/b[0])
 return F(out[n])
def arr0(a):return ar([coef(x,0)for x in a.ravel()]).reshape(a.shape)
def h_k(base,p,q):
 bj=killing(p,q,Bco);cj=killing(p,q,Cco);h=tuple(x+y for x,y in zip(bj,cj));phi=kbbjet(base,p,q);bc=symprod(base,bj,cj);k=tuple(x+Q(8,5)*y for x,y in zip(phi,bc));return h,k

def base_hessian(base):
 R=base.R;X=eye(6)[0];Y=eye(6)[3];N=eye(6)[[1,2,4,5]]
 HA=2*eins('ijkl,ai,j,k,bl->ab',R,N,Y,Y,N);HB=2*eins('ijkl,i,aj,bk,l->ab',R,X,N,N,X)
 HAB=2*(eins('ijkl,ai,bj,k,l->ab',R,N,N,Y,X)+eins('ijkl,ai,j,bk,l->ab',R,N,Y,N,X))
 H=np.block([[HA,HAB],[HAB.T,HB]])
 rg=np.r_[2*eins('ijkl,ai,j,k,l->a',R,N,Y,Y,X),2*eins('ijkl,i,aj,k,l->a',R,X,N,Y,X)]
 return H,rg

def extract_multilinear(comp,vs,n):
 # vs list of two-term formal vector polynomials [constant, linear].
 out=zero
 import itertools
 for degs in itertools.product(range(2),repeat=4):
  d=sum(degs)
  if d>n:continue
  for (i,x),(j,y),(k,zv),(l,w) in itertools.product(*[nonzero(vs[k][degs[k]])for k in range(4)]):
   out+=coef(comp(i,j,k,l),n-d)*x*y*zv*w
 return out
