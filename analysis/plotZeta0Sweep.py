#!/usr/bin/env python3
"""Plot whatever zeta_0 arms are on disk. Re-runnable: skips missing arms."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, numpy as np, glob, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matio          # reads -v7 and -v7.3 alike; scipy cannot open -v7.3
SURF='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; GRID='#e3e2df'
RAMP=['#9ec5f4','#5598e7','#2a78d6','#184f95','#0d366b']       # sequential, low->high zeta_0
VALS=[('0p1','0.1'),('1','1.0'),('10','10'),('100','100 (manuscript)'),('1000','1000')]
def load(tag):
    fs=sorted(glob.glob(f'analysis/probes/data/chain/chain_{tag}_leg*.mat'))
    if not fs: return None
    T=[];B=[];S=[];E=[];P=[]
    for i,f in enumerate(fs):
        d=matio.load_rec(f)
        g=lambda k: np.atleast_1d(np.asarray(d[k]).squeeze())
        ts=g('ts')
        # Files written before the 2026-08-25 fix double-counted the leg offset
        # (leg 2 stored 20..30 instead of 10..20). Detect rather than assume:
        # leg i (0-based) should start at i*10.
        if abs(ts[0]-(i+1)*10.0) < abs(ts[0]-i*10.0): ts = ts - 10.0
        k=slice(1,None) if i>0 else slice(None)
        T.append(ts[k]); B.append(g('bedInt')[k]); S.append(g('supInt')[k])
        E.append(g('effO2')[k]*1e3); P.append(np.asarray(d['phiT'])[:,k])
        z=g('z'); z0=float(np.ravel(d['zeta0'])[0])
    return dict(ts=np.concatenate(T),bed=np.concatenate(B),sup=np.concatenate(S),
                eff=np.concatenate(E),phi=np.concatenate(P,axis=1),z=z,z0=z0,legs=len(fs))
def run(N):
    arms=[(lab,c,load(f'z0_{t}_n{N}')) for (t,lab),c in zip(VALS,RAMP)]
    arms=[a for a in arms if a[2] is not None]
    if not arms: print(f'N={N}: nothing on disk yet'); return
    fig,ax=plt.subplots(1,4,figsize=(19,4.4),facecolor=SURF)
    for a in ax:
        a.set_facecolor(SURF)
        for s in ('top','right'): a.spines[s].set_visible(False)
        for s in ('left','bottom'): a.spines[s].set_color(GRID)
        a.grid(True,color=GRID,lw=0.8,zorder=0); a.set_axisbelow(True); a.tick_params(colors=INK2,labelsize=9)
    for lab,c,r in arms:
        ax[0].plot(r['ts'],r['bed'],color=c,lw=2,zorder=3,label=f'$\\zeta_0$ = {lab}')
        k=int(np.argmax(r['bed'])); ax[0].plot([r['ts'][k]],[r['bed'][k]],'o',color=c,ms=8,zorder=4)
        ax[1].plot(r['ts'],r['sup'],color=c,lw=2,zorder=3)
        ax[2].plot(r['phi'][:,-1],r['z'],color=c,lw=2,zorder=3)
        ax[3].plot(r['ts'],r['eff'],color=c,lw=2,zorder=3)
    ax[0].set_title('bed biomass',color=INK,fontsize=11,loc='left',pad=8)
    ax[0].set_ylabel('$\\int\\varepsilon\\phi_b dz$, z $\\geq-\\delta$  [m]',color=INK2,fontsize=10)
    ax[0].legend(frameon=False,fontsize=9,labelcolor=INK2)
    ax[1].set_title('supernatant biomass',color=INK,fontsize=11,loc='left',pad=8)
    ax[1].set_ylabel('$\\int\\varepsilon\\phi_b dz$, z $<-\\delta$  [m]',color=INK2,fontsize=10)
    ax[2].set_title('final profile',color=INK,fontsize=11,loc='left',pad=8)
    ax[2].set_xlabel('$\\phi_b$',color=INK2,fontsize=10); ax[2].set_ylabel('z [m]',color=INK2,fontsize=10)
    ax[2].set_ylim(0.6,-0.4); ax[2].axhline(0,color=INK2,lw=1,ls='--',zorder=2)
    ax[3].set_title('effluent O$_2$',color=INK,fontsize=11,loc='left',pad=8)
    ax[3].set_ylabel('effluent O$_2$  [g/m³]',color=INK2,fontsize=10)
    for a in (ax[0],ax[1],ax[3]): a.set_xlabel('t [d]',color=INK2,fontsize=10)
    tmax=max(r['ts'][-1] for _,_,r in arms)
    fig.suptitle(f'$\\zeta_0$ sweep — N = {N}, $\\zeta_1$ = 0.27, to t = {tmax:.0f} d '
                 f'({len(arms)}/5 arms on disk)',color=INK,fontsize=13,x=0.006,ha='left',y=0.98)
    fig.tight_layout(rect=[0,0,1,0.92])
    out=f'analysis/results/figures/zeta0_sweep_n{N}.png'
    fig.savefig(out,dpi=150,facecolor=SURF); print('wrote',out)
    print(f"  {'zeta_0':>8s}{'legs':>6s}{'t_end':>7s}{'bed peak':>10s}{'at t':>6s}{'decline':>9s}{'bed(end)':>10s}{'supFrac':>9s}{'effO2':>8s}")
    for lab,_,r in arms:
        k=int(np.argmax(r['bed'])); sf=100*r['sup'][-1]/(r['sup'][-1]+r['bed'][-1])
        print(f"  {lab:>8s}{r['legs']:6d}{r['ts'][-1]:7.1f}{r['bed'][k]:10.5f}{r['ts'][k]:6.1f}"
              f"{r['bed'][k]/max(r['bed'][-1],1e-30):8.2f}x{r['bed'][-1]:10.5f}{sf:8.1f}%{r['eff'][-1]:8.3f}")
for N in (sys.argv[1:] or ['200','500']): run(N)
