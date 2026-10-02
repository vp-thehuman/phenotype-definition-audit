#!/usr/bin/env python3
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
O=f"{ROOT}/mr/out"
INK="#1a1a1a"; MUTED="#6b6b6b"; GRID="#e3e3e3"; BLUE="#2f6f8f"; RED="#a8443a"; GREY="#9aa5ab"; NULLC="#b0b0b0"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.edgecolor":GRID,
   "axes.linewidth":0.8,"text.color":INK,"xtick.color":MUTED,"ytick.color":INK,"figure.dpi":300})
ORDER=["E1_EAGLE_criteria","E3_UKB_ICD","E5_BUDU_meta","E2_UKB_selfreport","E4_allergic_broad"]
SHORT={"E1_EAGLE_criteria":"Clinician criteria","E3_UKB_ICD":"Hospital ICD",
       "E5_BUDU_meta":"Pooled meta-analysis","E2_UKB_selfreport":"Hay fever/rhinitis/eczema",
       "E4_allergic_broad":"Allergic composite"}

# ============ FIG 3: the mechanism, MVMR dumbbell ============
mv=pd.read_csv(f"{O}/Table3_mvmr.csv"); mv=mv[mv.outcome=="Asthma"]
NAME={"Hospital ICD (UK Biobank)":"Hospital ICD","Pooled meta-analysis":"Pooled meta-analysis",
      "UKB hayfever/rhinitis/eczema":"Hay fever/rhinitis/eczema","Broad allergic composite":"Allergic composite"}
mv["short"]=mv.definition.map(NAME)
mv=mv.set_index("short").loc[["Hospital ICD","Pooled meta-analysis","Hay fever/rhinitis/eczema","Allergic composite"]].reset_index()
fig,ax=plt.subplots(figsize=(10.4,4.1))
y=np.arange(len(mv))[::-1]
for yy,(_,r) in zip(y,mv.iterrows()):
    ax.plot([r.uni_OR,r.mvmr_OR_definition],[yy,yy],color=GREY,lw=2.2,zorder=1,solid_capstyle="round")
    ax.plot(r.uni_OR,yy,"o",ms=9,color=GREY,mec="white",mew=1.5,zorder=2)
    col=RED if r.mvmr_p<0.05 else BLUE
    ax.plot(r.mvmr_OR_definition,yy,"o",ms=10,color=col,mec="white",mew=1.6,zorder=3)
    ax.hlines(yy,r.mvmr_lo,r.mvmr_hi,color=col,lw=2.0,alpha=.45,zorder=2)
    txt="absorbed by AD liability" if r.mvmr_p>=0.05 else f"retains independent effect (p = {r.mvmr_p:.0e})".replace("e-0","e-")
    ax.text(3.18,yy,txt,va="center",ha="right",fontsize=8,color=col,
            weight="bold" if r.mvmr_p<0.05 else "normal")
ax.axvline(1.0,color=NULLC,lw=1.0,ls=(0,(4,3)),zorder=0)
ax.set_yticks(y); ax.set_yticklabels(mv["short"],fontsize=9.2); ax.tick_params(axis="y",length=0,pad=6)
ax.set_xlim(0.85,3.22); ax.set_ylim(-0.6,len(mv)-0.4)
ax.xaxis.set_major_locator(FixedLocator([1.0,1.25,1.5,1.75,2.0,2.25]))
ax.set_xlabel("MR odds ratio, AD to asthma",fontsize=8.8,color=MUTED)
ax.grid(axis="x",color=GRID,lw=0.7); ax.set_axisbelow(True)
for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
ax.set_title("Conditioning on clinician-criteria AD liability absorbs strict definitions\nbut not loose ones",
             fontsize=12.4,color=INK,loc="left",pad=12,weight="bold")
ax.text(0,1.005,"",transform=ax.transAxes)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker="o",ls="",color=GREY,ms=8,label="univariable"),
                   Line2D([],[],marker="o",ls="",color=RED,ms=8,label="conditional on clinician-criteria AD liability")],
          loc="upper left",frameon=False,fontsize=8.4,ncol=2,bbox_to_anchor=(0.0,-0.14))
fig.savefig(f"{O}/Fig3_mvmr.png",bbox_inches="tight",facecolor="white")
fig.savefig(f"{O}/Fig3_mvmr.pdf",bbox_inches="tight",facecolor="white")

# ============ FIG 4: predictability, two panels ============
ix=pd.read_csv(f"{O}/Table5_contamination_index.csv")
S=pd.read_csv(f"{O}/Table6_outcome_ladder.csv")
fig,(a1,a2)=plt.subplots(1,2,figsize=(12.2,4.7),gridspec_kw=dict(wspace=0.30))
s=ix.dropna(subset=["mvmr_conditional_OR"])
a1.scatter(s.excess_rg,s.mvmr_conditional_OR,s=95,color=RED,edgecolor="white",linewidth=1.5,zorder=3)
z=np.polyfit(s.excess_rg,np.log(s.mvmr_conditional_OR),1); xs=np.linspace(-0.12,0.44,50)
a1.plot(xs,np.exp(np.polyval(z,xs)),color=MUTED,lw=1.2,ls=(0,(5,3)),zorder=1)
for _,r in s.iterrows():
    a1.annotate(r.label,(r.excess_rg,r.mvmr_conditional_OR),textcoords="offset points",
                xytext=(0,13),ha="center",fontsize=8,color=INK)
a1.axhline(1.0,color=NULLC,lw=1.0,ls=(0,(4,3)))
a1.set_xlabel("Excess genetic correlation with the contaminant\n(observed minus that predicted by AD liability alone)",fontsize=8.5,color=MUTED)
a1.set_ylabel("MVMR odds ratio conditional on\nclinician-criteria AD liability",fontsize=8.5,color=MUTED)
a1.set_xlim(-0.14,0.46); a1.set_ylim(0.92,2.30)
a1.set_title("A. The index predicts which definitions are contaminated",fontsize=10.4,loc="left",pad=10,weight="bold")
r1=np.corrcoef(s.excess_rg,np.log(s.mvmr_conditional_OR))[0,1]
a1.text(0.97,0.06,f"r = {r1:.2f}   (k = {len(s)})",transform=a1.transAxes,ha="right",fontsize=9.4,color=RED,weight="bold")

a2.scatter(S.rg_with_rhinitis,S.spread,s=95,color=BLUE,edgecolor="white",linewidth=1.5,zorder=3)
z2=np.polyfit(S.rg_with_rhinitis,np.log(S.spread),1); xs2=np.linspace(-0.02,1.05,50)
a2.plot(xs2,np.exp(np.polyval(z2,xs2)),color=MUTED,lw=1.2,ls=(0,(5,3)),zorder=1)
OFF={"Allergic rhinitis":(-6,10,"right"),"Asthma":(8,8,"left"),"COPD":(8,6,"left"),
     "Migraine":(0,-17,"center"),"Anxiety":(6,8,"left"),"Depression":(-6,-16,"right"),"IBD":(8,4,"left")}
for _,r in S.iterrows():
    dx,dy,ha=OFF.get(r.outcome,(0,12,"center"))
    a2.annotate(r.outcome,(r.rg_with_rhinitis,r.spread),textcoords="offset points",
                xytext=(dx,dy),ha=ha,fontsize=8,color=INK)
a2.axhline(1.0,color=NULLC,lw=1.0,ls=(0,(4,3)))
a2.set_xlabel("Genetic correlation of the outcome\nwith the contaminating condition",fontsize=8.5,color=MUTED)
a2.set_ylabel("Spread of MR estimates across\nthe five AD definitions (max / min)",fontsize=8.5,color=MUTED)
a2.set_xlim(-0.06,1.12); a2.set_ylim(0.96,1.92)
a2.set_title("B. And predicts which outcomes are affected",fontsize=10.4,loc="left",pad=10,weight="bold")
r2=np.corrcoef(S.rg_with_rhinitis,np.log(S.spread))[0,1]
a2.text(0.97,0.06,f"r = {r2:.2f}   (k = {len(S)})",transform=a2.transAxes,ha="right",fontsize=9.4,color=BLUE,weight="bold")
for ax in (a1,a2):
    ax.grid(color=GRID,lw=0.7); ax.set_axisbelow(True)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
fig.suptitle("Definitional bias in MR is predictable from summary statistics alone",
             fontsize=12.8,color=INK,x=0.005,ha="left",y=1.03,weight="bold")
fig.savefig(f"{O}/Fig4_predictability.png",bbox_inches="tight",facecolor="white")
fig.savefig(f"{O}/Fig4_predictability.pdf",bbox_inches="tight",facecolor="white")
print("figures written")
