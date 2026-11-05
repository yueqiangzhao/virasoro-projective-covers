from contextlib import redirect_stdout
from io import StringIO
with redirect_stdout(StringIO()):
 from mixed_candidate import *
def P(vertices,cs):return tuple(idx[(vertices[k],vertices[k+1],cs[k])] for k in range(2))
rr=[]
def add(*terms):
 r={}
 for c,vertices,cs in terms:
  p=P(vertices,cs);r[p]=r.get(p,F(0))+c
 rr.append({p:c for p,c in r.items() if c})
for i,j in [('A','C'),('A','D'),('B','C'),('B','D'),('C','A'),('C','B'),('D','A'),('D','B')]:
 add((1,i+j+i,(1,1)));add((1,i+j+i,(-1,-1)))
 add((1,i+j+i,(1,-1)),(1,i+j+i,(-1,1)))
add((1,'TAT',(0,0)),(1,'TBT',(0,0)))
for i in ['A','B']:
 add((1,i+'T'+i,(0,0)),(1,i+'C'+i,(1,-1)),(1,i+'D'+i,(1,-1)))
for s,m1,m2,t in [('A','D','C','B'),('B','D','C','A'),('C','B','A','D'),('D','B','A','C')]:
 for u,v in product([1,-1],repeat=2):
  add((1,s+m1+t,(u,v)),(-1,s+m2+t,(v,u)))
for s,t in [('A','B'),('B','A')]:
 add((1,s+'T'+t,(0,0)),(-1,s+'C'+t,(1,-1)),(1,s+'C'+t,(-1,1)))
for u in [1,-1]:
 add((1,'TAC',(0,u)),(-1,'TBC',(0,u)))
 add((1,'CAT',(u,0)),(-1,'CBT',(u,0)))
 add((1,'TAD',(0,u)),(1,'TBD',(0,u)))
 add((1,'DAT',(u,0)),(1,'DBT',(u,0)))
pp=[p for s in V for p in paths(s,2)];ids={p:j for j,p in enumerate(pp)}
rank=len(rref([{ids[p]:c for p,c in r.items()} for r in rr]))
checks=0
for r in rr:
 p=next(iter(r));s=a[p[0]][0];t=end(p)
 for q in paths(t,2):
  if end(q)==s:
   checks+=1
   assert sum(c*phi(p+q) for p,c in r.items())==0,(r,q)
print('explicit relation vectors',len(rr),'rank',rank,'potential annihilator checks',checks)
assert len(rr)==rank==53
# N=sum N0,N1,N2, from closed formula in proof; squared values.
Ns={
'T':{P('TAT',(0,0)):F(1)},
'A':{P('ACA',(1,-1)):F(1),P('ADA',(1,-1)):F(-3)},
'B':{P('BCB',(1,-1)):F(1),P('BDB',(1,-1)):F(1)},
'C':{P('CAC',(1,-1)):F(1),P('CBC',(1,-1)):F(-3)},
'D':{P('DAD',(1,-1)):F(1),P('DBD',(1,-1)):F(1)}}
for ar,(s,t,_) in enumerate(a):
 for rev in out[t]:
  if a[rev][1]!=s:continue
  assert sum(c*phi(p+(ar,rev)) for p,c in Ns[s].items())==sum(c*phi((ar,)+p+(rev,)) for p,c in Ns[t].items())
print('N centrality checked against all reverse arrows')
print('N^2 socle trace coefficients',{s:str(sum(c*d*phi(p+q) for p,c in n.items() for q,d in n.items())) for s,n in Ns.items()})
