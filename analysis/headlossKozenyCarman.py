"""Kozeny-Carman head-loss diagnostic, post-processed from chained probe frames.

Reads chain_<tag>_leg*.mat (analysis/probes/data/chain), takes the pore fraction
eps_eff = eps*(1 - phi_b) on the uniform bed (z >= delta, where eps = 0.4) and
reports H/H0 = (1/L) int K0/K dz with K ~ eps^3/(1-eps)^2 -- the bed head loss
relative to the clean bed. The roughness ramp 0 <= z < delta (eps -> 1) has no
clean-bed reference and is excluded. Proposed in .claude/CRITIQUE.md section 7.1;
literature band for operational clogging is H/H0 = 2-10 (Kim2010, Volk2016).
"""
import sys; sys.path.insert(0,'analysis')
import matio, numpy as np, glob
def chain(tag):
    fs=sorted(glob.glob(f'analysis/probes/data/chain/chain_{tag}_leg*.mat'),key=lambda f:int(f.rsplit('leg',1)[1].split('.')[0]))
    T=[];PR=[];B=[];S=[]
    for i,f in enumerate(fs):
        d=matio.load_rec(f); g=lambda k: np.atleast_1d(np.asarray(d[k]).squeeze())
        ts=g('ts')
        if i>0 and abs(ts[0]-2*T[-1][-1])<abs(ts[0]-T[-1][-1]): ts=ts-T[-1][-1]
        k=slice(1,None) if i>0 else slice(None)
        T.append(ts[k]);PR.append(np.asarray(d['phiT'])[:,k]);B.append(g('bedInt')[k]);S.append(g('supInt')[k])
        z=np.ravel(d['z']); eps=np.ravel(d['eps']); delta=float(np.ravel(d['delta'])[0])
    return dict(ts=np.concatenate(T),phi=np.concatenate(PR,axis=1),bed=np.concatenate(B),sup=np.concatenate(S),z=z,eps=eps,delta=delta)
def headloss(r):
    # Kozeny-Carman on the uniform bed z >= delta (eps = 0.4). The roughness ramp
    # (0 <= z < delta, eps -> 1) has no clean-bed reference and is excluded.
    z=r['z']; sel=z>=r['delta']; e0=r['eps'][sel][:,None]
    e=np.clip(e0*(1-r['phi'][sel,:]),1e-9,None)
    K0_over_K=(e0/e)**3*((1-e)/(1-e0))**2
    dz=np.gradient(z)[sel][:,None]
    H=np.sum(K0_over_K*dz,axis=0)/np.sum(dz)      # total bed head loss / clean-bed head loss
    return H,K0_over_K
if __name__=='__main__':
    tags={1:'z0_1_n500',3:'z0w_3_n500',5:'z0w_5_n500_local',10:'z0w_10_n500',30:'z0w_30_n500',100:'z0_100_n500'}
    print("Kozeny-Carman: H/H0 = bed head loss relative to clean bed (z >= delta only)")
    print(f"{'z0':>4s}{'t_end':>6s}{'H/H0 5d':>9s}{'10d':>8s}{'20d':>8s}{'end':>9s}{'t(2x)':>8s}{'t(3x)':>8s}{'t(10x)':>8s}{'top-cell K0/K end':>19s}")
    for z0,tag in sorted(tags.items()):
        r=chain(tag); H,R=headloss(r); ts=r['ts']
        def at(t):
            j=int(np.argmin(abs(ts-t))); return f"{H[j]:8.2f}" if abs(ts[j]-t)<0.3 else f"{'--':>8s}"
        def tc(x):
            i=int(np.argmax(H>=x)); return f"{ts[i]:8.1f}" if H[i]>=x else f"{'>end':>8s}"
        print(f"{z0:4d}{ts[-1]:6.0f}{at(5)}{at(10)}{at(20)}{H[-1]:9.2f}{tc(2)}{tc(3)}{tc(10)}{R[0,-1]:19.1f}")
