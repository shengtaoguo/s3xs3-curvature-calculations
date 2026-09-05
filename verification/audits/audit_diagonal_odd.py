"""Full base-and-plane audit of the odd quadratic term at the diagonal."""
from exact_engine import *
from exact_sparse import Linear
from exact_diagonal import coef,arr0,base_hessian
p=eye(4)[0];X=eye(6)[0];Y=eye(6)[3];N=eye(6)[[1,2,4,5]]
H10=zeros((10,10));r10=zeros(10)
for ax in (2,3):
    q=ar([(1-z*z)/(1+z*z),zero,zero,zero]);q[ax]=2*z/(1+z*z)
    base=ExactBase(*g_jets(q));base.X=X;base.Y=Y;base.N=N
    cs=coord_jets(p,q)
    hj=tuple(a+b for a,b in zip(const_height(cs,2,AP),const_height(cs,7,AQ)))
    lin=Linear(base,hj);ss,rr=lin.data();H8,rg=base_hessian(base)
    H10[:8,:8]=arr0(H8);idx=ax+6
    H10[:8,idx]=H10[idx,:8]=ar([coef(x,1)for x in rg]);H10[idx,idx]=2*coef(base.R[0,3,3,0],2)
    r10[:8]=arr0(rr);r10[idx]=coef(lin.component(0,3,3,0),1)
inds=[0,1,2,3,4,5,8,9];hh=H10[np.ix_(inds,inds)];rh=r10[inds]
eta=zeros(10);eta[inds]=-inv(hh)@rh
assert all(x==0 for x in H10@eta+r10)
g=arr0(base.g);s=arr0(ss);D=g[0,0]*g[3,3]-g[0,3]**2
fixed=s[1]@g@s[1]-s[0]@g@s[2]
value=(fixed-rh@inv(hh)@rh/2)/D
expected=Q(12632468,31752000*25)
assert value==expected
print('PASS: full base-and-plane odd quadratic term at the diagonal:',value)
