from exact_engine import *
from exact_sparse import Linear
from exact_cubic import symprod
pairs6=[(i,j)for i in range(6)for j in range(i,6)]
# Coefficients of the explicit tensor K_BO.  These are exact rational numbers;
# no numerical fit or external data file is used in this certificate.
boco=zeros((16,21))
def add(group,pair,value):
 boco[group,pairs6.index(tuple(sorted(pair)))]+=value

def A_group(group,factor):
 for pair,co in [((0,1),one),((0,4),-one),((1,3),one),((3,4),-one)]:add(group,pair,factor*co)
def W_group(group,factor,sign):
 for pair,co in [((0,2),Q(13,7)),((0,5),Q(sign)),((2,3),Q(13,7)),((3,5),Q(sign))]:add(group,pair,factor*co)
def T_group(group,factor):
 for pair in [(0,0),(0,3),(3,3)]:add(group,pair,-Q(4,7)*factor)
a=Q(11,325);b=Q(6,65)
A_group(0,a);W_group(1,a,-1);T_group(3,a)
W_group(4,a,1)
for pair in [(0,1),(0,4),(1,3),(3,4)]:add(5,pair,-a)
T_group(6,a)
W_group(8,b,1)
for pair in [(0,1),(0,4),(1,3),(3,4)]:add(9,pair,-b)
T_group(10,b)
A_group(12,b);W_group(13,b,-1);T_group(15,b)

def rotate_deriv(p,pu,a,b):
 A=zeros((4,4));A[a,b]=-one;A[b,a]=one
 dp=A@p;dpu=A@pu;pc=np.r_[p[0],-p[1:]];dpc=np.r_[dp[0],-dp[1:]]
 u=qmul(pc,pu)[1:];w=qmul(pu,pc)[1:]
 du=(qmul(dpc,pu)+qmul(pc,dpu))[1:];dw=(qmul(dpu,pc)+qmul(pu,dpc))[1:]
 return A,dp,du,dw,u,w

def bo_Lgrad(p,pu,q,r0,a,b,D):
 A,dp,du,dw,u,w=rotate_deriv(p,pu,a,b);dq=A@q
 zz=np.r_[u,w];dz=np.r_[du,dw]
 zzp=ar([zz[i]*zz[j]for i,j in pairs6]);dzp=ar([dz[i]*zz[j]+zz[i]*dz[j]for i,j in pairs6])
 scal=np.r_[p,q,r0*p,r0*q];ds=np.r_[dp,dq,r0*dp,r0*dq]
 val=eins('a,b,ab->',scal,zzp,boco).item()/(2*D)
 deriv=(eins('a,b,ab->',ds,zzp,boco).item()+eins('a,b,ab->',scal,dzp,boco).item())/(2*D)
 return val,deriv

def L_on_flat(jet,base):
 _,_,dd=jet;X,Y=base.X,base.Y
 return (eins('abij,a,b,i,j->',dd,X,Y,X,Y).item()-eins('abij,a,b,i,j->',dd,X,X,Y,Y).item()/2-eins('abij,a,b,i,j->',dd,Y,Y,X,X).item()/2)/base.D
