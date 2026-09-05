from exact_engine import *
from exact_sparse import Linear
from exact_cubic import scalar_prod,scalar_times

def relative_jets(r):
 dr=zeros((6,4));ddr=zeros((6,6,4))
 for i in range(3):
  dr[i]=-Lquat[i]@r;dr[i+3]=Rquat[i]@r
  for l in range(3):
   ddr[l,i]=Lquat[i]@Lquat[l]@r;ddr[l+3,i+3]=Rquat[i]@Rquat[l]@r
   ddr[l,i+3]=-Rquat[i]@Lquat[l]@r;ddr[l+3,i]=-Lquat[i]@Rquat[l]@r
 return r,dr,ddr

def h0jet(r):
 cs=relative_jets(r);P=zeros((3,3));dP=zeros((6,3,3));ddP=zeros((6,6,3,3));Pj=[P,dP,ddP]
 for i in range(3):
  for j in range(3):
   sc=scalar_prod(tuple(x[...,i+1]for x in cs),tuple(x[...,j+1]for x in cs))
   for a in range(3):Pj[a][...,i,j]-=sc[a]
   if i==j:
    for k in range(3):
     for a in range(3):Pj[a][...,k,k]+=sc[a]
 mats=[zeros((6,6)),zeros((6,6,6)),zeros((6,6,6,6))]
 for i in range(3):mats[i][...,:3,3:]=Pj[i];mats[i][...,3:,:3]=Pj[i]
 return scalar_times(tuple(x[...,0]for x in cs),tuple(mats))

def at0(x):
 if isinstance(x,FracElement):return F(x.numer.get((0,),QQ.zero)/x.denom.get((0,),QQ.zero))
 return F(x)
def deriv0(x,n=1):
 for i in range(n):x=x.diff(z)
 return at0(x)
