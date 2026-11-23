import numpy as np
exec(open("t3_1d.py").read().split("for nm,M in")[0])
def run_seq(M, dts):
    Tn=np.zeros(N); t=0.0
    Mff,Kff=M[np.ix_(fr,fr)],K[np.ix_(fr,fr)]; MfB,KfB=M[fr,N-1],K[fr,N-1]
    for dt in dts:
        b0,b1=gB(t),gB(t+dt); t+=dt
        A=Mff/dt+Kff; rhs=Mff@Tn[fr]/dt-MfB*(b1-b0)/dt-KfB*b1
        Tnew=Tn.copy(); Tnew[fr]=np.linalg.solve(A,rhs); Tnew[N-1]=b1; Tn=Tnew
    return Tn[8], t
def seq(first, growth, cap, T=32.0):
    out=[]; t=0.0; dt=first
    while t < T-1e-12:
        d=min(dt, T-t); out.append(d); t+=d; dt=min(dt*growth, cap)
    return out
from scipy.optimize import brentq
f=lambda dt: run_seq(Mc,[dt]*int(round(32/dt)) if abs(32/dt-round(32/dt))<1e-9 else seq(dt,1,dt))[0]-37.25446
for dt in (0.1,0.2,0.25,0.32,0.4,0.5): print("const dt",dt, round(run_seq(Mc, seq(dt,1.0,dt))[0],4))
for g in (1.25,1.5,2.0):
    s=seq(0.1,g,0.5); print("start 0.1 growth",g,"cap 0.5: n_inc",len(s), round(run_seq(Mc,s)[0],4))
