import numpy as np
L,n,k,rc=0.1,10,35.0,7200*440.5
h=L/n; N=n+1
K=np.zeros((N,N)); Mc=np.zeros((N,N))
for e in range(n):
    i=[e,e+1]
    K[np.ix_(i,i)]+=k/h*np.array([[1,-1],[-1,1]])
    Mc[np.ix_(i,i)]+=rc*h/6*np.array([[2,1],[1,2]])
Ml=np.diag(Mc.sum(1))
gB=lambda t:100*np.sin(np.pi*t/40)
fr=np.arange(1,N-1)  # free nodes (0 and N-1 Dirichlet)
def run(M,scheme,dt,T=32.0):
    nst=int(round(T/dt)); Tn=np.zeros(N); Tp=None
    Mff,Kff=M[np.ix_(fr,fr)],K[np.ix_(fr,fr)]
    MfB,KfB=M[fr,N-1],K[fr,N-1]
    for s in range(nst):
        t1=(s+1)*dt; b1=gB(t1); b0=gB(s*dt)
        if scheme=='BE' or (scheme=='BDF2' and s==0):
            A=Mff/dt+Kff
            rhs=Mff@Tn[fr]/dt - MfB*(b1-b0)/dt - KfB*b1
            Tnew=Tn.copy(); Tnew[fr]=np.linalg.solve(A,rhs); Tnew[N-1]=b1
        else: # BDF2
            bm=gB((s-1)*dt)
            A=1.5*Mff/dt+Kff
            rhs=Mff@(2*Tn[fr]-0.5*Tp[fr])/dt - MfB*(1.5*b1-2*b0+0.5*bm)/dt - KfB*b1
            Tnew=Tn.copy(); Tnew[fr]=np.linalg.solve(A,rhs); Tnew[N-1]=b1
        Tp,Tn=Tn,Tnew
    return Tn[8]
for nm,M in [('consistent',Mc),('lumped',Ml),('mixed 1/2(Mc+Ml)',0.5*(Mc+Ml))]:
    r=[run(M,'BDF2',0.0125)]+[run(M,'BE',d) for d in (0.1,0.25,0.5)]
    print(f"{nm:18s} BDF2(0.0125)={r[0]:.4f} ({(r[0]/36.6-1)*100:+.2f}%)  BE0.1={r[1]:.4f} ({(r[1]/36.6-1)*100:+.2f}%) BE0.25={r[2]:.4f} BE0.5={r[3]:.4f} ({(r[3]/36.6-1)*100:+.2f}%)")
