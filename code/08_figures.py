#!/usr/bin/env python3
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator

r=pd.read_csv(f"{ROOT}/mr/out/mr_results_v2.csv"); het=pd.read_csv(f"{ROOT}/mr/out/heterogeneity_v2.csv")
con=pd.read_csv(f"{ROOT}/mr/out/Table2_definition_concordance.csv")
ORDER=["E1_EAGLE_criteria","E3_UKB_ICD","E5_BUDU_meta","E2_UKB_selfreport","E4_allergic_broad"]
LAB={"E1_EAGLE_criteria":"Clinician criteria\nEAGLE, 18,900 cases",
     "E3_UKB_ICD":"Hospital ICD record\nUK Biobank, 12,176 cases",
     "E5_BUDU_meta":"Pooled meta-analysis\n2023, 60,653 cases",
     "E2_UKB_selfreport":"Hay fever, rhinitis or eczema\nUK Biobank touchscreen, n = 461,199",
     "E4_allergic_broad":"Broad allergic composite\n180,129 cases"}
SHORT={"E1_EAGLE_criteria":"Clinician criteria","E3_UKB_ICD":"Hospital ICD","E5_BUDU_meta":"Pooled meta-analysis",
       "E2_UKB_selfreport":"Hay fever/rhinitis/eczema","E4_allergic_broad":"Allergic composite"}
INK="#1a1a1a"; MUTED="#6b6b6b"; GRID="#e3e3e3"; BLUE="#2f6f8f"; RED="#a8443a"; NULLC="#b0b0b0"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.edgecolor":GRID,
                     "axes.linewidth":0.8,"text.color":INK,"xtick.color":MUTED,"ytick.color":INK,"figure.dpi":300})

# ---------------- Figure 1: forest, mental-health outcomes + asthma positive control
PAN=[("O1_FG_depression","Depression\nFinnGen R11",BLUE,(0.90,1.245),[0.95,1.00,1.05,1.10]),
     ("O2_FG_anxiety","Anxiety\nFinnGen R11",BLUE,(0.90,1.245),[0.95,1.00,1.05,1.10]),
     ("O3_PGC_MDD","MDD, sensitivity\nPGC 2018 ex-23andMe",BLUE,(0.90,1.245),[0.95,1.00,1.05,1.10]),
     ("O4_FG_asthma","Asthma\nPOSITIVE CONTROL",RED,(1.02,2.30),[1.2,1.4,1.6,1.8])]
fig,axes=plt.subplots(1,4,figsize=(16.4,4.7),sharey=True,
    gridspec_kw=dict(wspace=0.11,left=0.155,right=0.99,top=0.745,bottom=0.185))
y=np.arange(len(ORDER))[::-1]
for ax,(ok,olab,col,xlim,ticks) in zip(axes,PAN):
    s=r[r.outcome==ok].set_index("exposure").reindex(ORDER)
    ax.axvline(1.0,color=NULLC,lw=1.0,ls=(0,(4,3)),zorder=1)
    ax.hlines(y,s.ivw_lo,s.ivw_hi,color=col,lw=2.0,alpha=.55,zorder=2)
    ax.plot(s.ivw_or,y,"o",ms=8,color=col,mec="white",mew=1.6,zorder=3)
    for yy,o_,lo,hi in zip(y,s.ivw_or,s.ivw_lo,s.ivw_hi):
        ax.text(xlim[1]-0.008*(xlim[1]-xlim[0]),yy,f"{o_:.2f} ({lo:.2f}-{hi:.2f})",
                va="center",ha="right",fontsize=7.5,color=MUTED,zorder=4)
    ax.set_xlim(*xlim); ax.set_ylim(-0.7,len(ORDER)-0.3)
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.set_xlabel("OR per unit higher log-odds of AD",fontsize=8.1,color=MUTED)
    ax.set_title(olab,fontsize=9.5,color=col if ok=="O4_FG_asthma" else INK,pad=8,
                 weight="bold" if ok=="O4_FG_asthma" else "normal")
    ax.grid(axis="x",color=GRID,lw=0.7,zorder=0); ax.set_axisbelow(True)
    for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
    h=het[het.outcome==ok].iloc[0]
    ax.text(0.42,-0.255,f"$I^2$ = {h.I2:.0f}%   p = {h.Q_p:.1e}   spread {h.max_min_ratio:.2f}x" if h.I2>20
            else f"$I^2$ = {h.I2:.0f}%   p = {h.Q_p:.2f}   spread {h.max_min_ratio:.2f}x",
            transform=ax.transAxes,ha="center",fontsize=7.8,color=col if h.I2>20 else MUTED,
            weight="bold" if h.I2>20 else "normal")
axes[0].set_yticks(y); axes[0].set_yticklabels([LAB[e] for e in ORDER],fontsize=8.1,linespacing=1.35)
axes[0].tick_params(axis="y",length=0,pad=6)
for a in axes[1:]: a.tick_params(axis="y",length=0,labelleft=False)
fig.suptitle("How atopic dermatitis is defined changes the Mendelian randomisation estimate,\n"
             "but only where there is a causal effect to distort",
             fontsize=13,color=INK,x=0.006,ha="left",y=1.005,weight="bold")
fig.text(0.006,0.845,"Same outcome, same pipeline, identical clumping (r²<0.001, 10 Mb, 1000G EUR). Only the exposure case definition varies. "
         "Null across mental-health outcomes; strongly definition-dependent for asthma.",fontsize=8.3,color=MUTED,ha="left")
fig.savefig(f"{ROOT}/mr/out/Fig1_forest.png",bbox_inches="tight",facecolor="white")
fig.savefig(f"{ROOT}/mr/out/Fig1_forest.pdf",bbox_inches="tight",facecolor="white")

# ---------------- Figure 2: liability capture vs estimate inflation
a=r[r.outcome=="O4_FG_asthma"].merge(con,left_on="exposure",right_on="definition")
fig2,ax=plt.subplots(figsize=(7.0,5.0))
ax.errorbar(a.lambda_vs_EAGLE,a.ivw_or,yerr=[a.ivw_or-a.ivw_lo,a.ivw_hi-a.ivw_or],
            fmt="o",ms=9,color=RED,ecolor=RED,elinewidth=1.6,alpha=.85,mec="white",mew=1.5,capsize=0)
z=np.polyfit(a.lambda_vs_EAGLE,np.log(a.ivw_or),1); xs=np.linspace(0.24,1.06,50)
ax.plot(xs,np.exp(np.polyval(z,xs)),color=MUTED,lw=1.2,ls=(0,(5,3)),zorder=0)
OFF={"E4_allergic_broad":(10,12,"left"),"E2_UKB_selfreport":(10,12,"left"),
     "E5_BUDU_meta":(0,14,"center"),"E1_EAGLE_criteria":(-8,16,"right"),
     "E3_UKB_ICD":(-10,-18,"right")}
for _,x_ in a.iterrows():
    dx,dy,ha=OFF[x_.exposure]
    ax.annotate(SHORT[x_.exposure],(x_.lambda_vs_EAGLE,x_.ivw_or),textcoords="offset points",
                xytext=(dx,dy),ha=ha,fontsize=8.1,color=INK)
ax.axhline(1.0,color=NULLC,lw=1.0,ls=(0,(4,3)))
ax.set_xlabel("$\\lambda$: strength of liability capture relative to clinician criteria\n(Deming slope of instrument effects vs EAGLE)",
              fontsize=8.6,color=MUTED)
ax.set_ylabel("MR odds ratio, AD → asthma",fontsize=8.6,color=MUTED)
ax.set_xlim(0.18,1.16); ax.set_ylim(0.98,2.02)
ax.grid(color=GRID,lw=0.7); ax.set_axisbelow(True)
for sp in ("top","right"): ax.spines[sp].set_visible(False)
rr=np.corrcoef(a.lambda_vs_EAGLE,np.log(a.ivw_or))[0,1]
ax.set_title("Looser definitions capture liability more weakly\nand inflate the MR estimate",
             fontsize=11.6,color=INK,loc="left",pad=10,weight="bold")
ax.text(0.98,0.05,f"r = {rr:.2f}",transform=ax.transAxes,ha="right",fontsize=10,color=RED,weight="bold")
fig2.savefig(f"{ROOT}/mr/out/Fig2_lambda_vs_estimate.png",bbox_inches="tight",facecolor="white")
fig2.savefig(f"{ROOT}/mr/out/Fig2_lambda_vs_estimate.pdf",bbox_inches="tight",facecolor="white")
print("figures written")
