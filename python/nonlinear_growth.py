# Direct nonlinear pseudo-spectral run of the four-field Hall system (Lockwood eqs 4.3-4.6) around psi0=cos x.
# Usage: python3 nonlinear_growth.py k eta h T   (e.g. 0.5 0.05 0.5 70 gives growth 0.1437; 1.2 0.05 0.5 120 decays)
import numpy as np, sys
def run(k,eta,h,T,Nx=64,Ny=32,dt=0.01,amp=1e-7,seed=1):
    Lx,Ly=2*np.pi,2*np.pi/k
    x=np.arange(Nx)*Lx/Nx; y=np.arange(Ny)*Ly/Ny; X,Y=np.meshgrid(x,y,indexing='ij')
    kx=np.fft.fftfreq(Nx,1/Nx)*2*np.pi/Lx; ky=np.fft.rfftfreq(Ny,1/Ny)*2*np.pi/Ly
    KX,KY=np.meshgrid(kx,ky,indexing='ij'); K2=KX**2+KY**2; K2i=np.where(K2>0,1/np.where(K2>0,K2,1),0)
    deal=(np.abs(KX)<(2/3)*np.abs(kx).max()+1e-9)&(np.abs(KY)<(2/3)*ky.max()+1e-9)
    F=lambda f: np.fft.rfft2(f); I=lambda g: np.fft.irfft2(g,s=(Nx,Ny))
    def br(fh,gh):
        return F(I(1j*KX*fh)*I(1j*KY*gh)-I(1j*KY*fh)*I(1j*KX*gh))*deal
    Spsi=F(eta*np.cos(X))
    rng=np.random.default_rng(seed)
    def noise():
        return amp*sum(rng.standard_normal()*np.cos(m*X+k*Y+rng.uniform(0,6.3)) for m in range(4))
    psi=F(np.cos(X)+noise()); om=F(noise()); w=F(noise()); b=F(noise())
    def rhs(s):
        om,psi,w,b=s; phi=-om*K2i; J=-K2*psi
        return (-br(phi,om)+br(psi,J),
                br(psi,phi+h*b)+eta*J+Spsi,
                -br(phi,w)+br(psi,b),
                -br(phi,b)+br(psi,w)-h*br(psi,J)-eta*K2*b)
    s=[om,psi,w,b]; iy=1  # ky=k column
    ts=[];amps=[]
    n=int(round(T/dt))
    for i in range(n+1):
        if i%int(round(1/dt))==0:
            e=sum(np.sum(np.abs(f[:,iy])**2) for f in s); ts.append(i*dt); amps.append(np.sqrt(e))
        k1=rhs(s); k2=rhs([a+dt/2*c for a,c in zip(s,k1)]); k3=rhs([a+dt/2*c for a,c in zip(s,k2)]); k4=rhs([a+dt*c for a,c in zip(s,k3)])
        s=[a+dt/6*(p+2*q+2*r+u) for a,p,q,r,u in zip(s,k1,k2,k3,k4)]
    ts=np.array(ts); la=np.log(np.array(amps))
    return ts,la
if __name__=="__main__":
    k,eta,h,T=map(float,sys.argv[1:5])
    ts,la=run(k,eta,h,T)
    m=len(ts); a=m//2
    g=np.polyfit(ts[a:],la[a:],1)[0]
    print(f"k={k} eta={eta} h={h}: fitted growth over t in [{ts[a]:.0f},{ts[-1]:.0f}] = {g:.5f}; amplitude {np.exp(la[0]):.2e} -> {np.exp(la[-1]):.2e}")
