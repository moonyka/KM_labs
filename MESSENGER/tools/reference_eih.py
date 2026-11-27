# Независимый покомпонентный эталон (4-26), без вызова prepareEIH.
import numpy as np
from pathlib import Path
positions=np.array([[0,0,0],[3,4,0],[-6,5,4]],dtype=np.longdouble)
velocities=np.array([[1,2,3],[-1,1,2],[3,-2,1]],dtype=np.longdouble)
mus=np.array([12,20,30],dtype=np.longdouble)
def ref(ri,vi,excluded=-1):
 ri=np.array(ri,dtype=np.longdouble);vi=np.array(vi,dtype=np.longdouble)
 def norm(r):return np.sqrt(sum(x*x for x in r))
 U_i=sum(mus[j]/norm(ri-positions[j]) for j in range(3) if j!=excluded)
 total=np.zeros(3,dtype=np.longdouble)
 for j in range(3):
  if j==excluded:continue
  aj=sum((mus[k]*(positions[k]-positions[j])/norm(positions[k]-positions[j])**3 for k in range(3) if k!=j),start=np.zeros(3,dtype=np.longdouble))
  U_j=sum(mus[k]/norm(positions[j]-positions[k]) for k in range(3) if k!=j)
  dij=ri-positions[j];dist=norm(dij);vj=velocities[j]
  q=-4*U_i-U_j+sum(vi*vi)+2*sum(vj*vj)-4*sum(vi*vj)-np.longdouble(1.5)*(sum(dij*vj)/dist)**2-sum(dij*aj)/2
  for c in range(3):
   total[c]+=mus[j]*(-dij[c])/dist**3*q
   total[c]+=mus[j]/dist**3*sum(dij*(4*vi-3*vj))*(vi[c]-vj[c])
   total[c]+=np.longdouble(3.5)*mus[j]*aj[c]/dist
 return total/np.longdouble('299792.458')**2
for r,v,ex in [([8,9,10],[2,3,4],-1),(positions[0],velocities[0],0),(positions[1],velocities[1],1)]:print(', '.join(format(float(x),'.18e') for x in ref(r,v,ex)))
