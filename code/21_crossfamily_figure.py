#!/usr/bin/env python3
import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
O=f"{ROOT}/mr/out"
INK="#1a1a1a"; MUTED="#6b6b6b"; GRID="#e3e3e3"; BLUE="#2f6f8f"; RED="#a8443a"; NULLC="#b0b0b0"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.edgecolor":GRID,
 "axes.linewidth":0.8,"text.color":INK,"xtick.color":MUTED,"ytick.color":INK,"figure.dpi":300})
ad=pd.read_csv(f"{O}/Table5_contamination_index.csv")
dep=pd.read_csv(f"{O}/Table7_depression_family_index.csv")
rc_dep=float(dep.rg_cidi_neuroticism.iloc[0])
rc_ad=float(ad[ad.trait=="AD_criteria"].rg_contaminant.iloc[0])

fig,(a1,a2)=plt.subplots(1,2,figsize=(12.4,5.0),gridspec_kw=dict(wspace=0.26,top=0.90))
xs=np.linspace(0,1.05,50)

# --- AD family
a1.plot(xs,xs*rc_ad,color=MUTED,lw=1.3,ls=(0,(5,3)),zorder=1)
a1.fill_between(xs,xs*rc_ad,1.05,color=RED,alpha=0.055,zorder=0)
a1.scatter(ad.rg_vs_criteria,ad.rg_contaminant,s=110,color=RED,edgecolor="white",lw=1.5,zorder=3)
OFF={"Clinician criteria":(-11,-14,"right"),"Hospital ICD":(11,-3,"left"),
     "Pooled meta-analysis":(-11,6,"right"),"Hay fever/rhinitis/eczema":(8,-12,"left"),"Allergic composite":(-10,7,"right")}
for _,r in ad.iterrows():
    dx,dy,ha=OFF.get(r.label,(0,12,"center"))
    a1.annotate(r.label,(r.rg_vs_criteria,r.rg_contaminant),textcoords="offset points",
                xytext=(dx,dy),ha=ha,fontsize=8,color=INK)
a1.set_title("A. Atopic dermatitis\ncontamination tracks definition breadth",fontsize=10.6,loc="left",pad=10,weight="bold")
a1.set_xlabel("Genetic correlation with the reference definition\n(clinician criteria)",fontsize=8.6,color=MUTED)
a1.set_ylabel("Genetic correlation with the contaminant\n(allergic rhinitis)",fontsize=8.6,color=MUTED)
a1.set_xlim(0.52,1.10); a1.set_ylim(0.28,1.02)

# --- Depression family
d=dep.copy()
d=pd.concat([d,pd.DataFrame([dict(label="CIDI lifetime depression\n(reference)",rg_vs_reference=1.0,
                                  rg_contaminant=rc_dep,excess_rg=0.0)])],ignore_index=True)
a2.plot(xs,xs*rc_dep,color=MUTED,lw=1.3,ls=(0,(5,3)),zorder=1)
a2.fill_between(xs,xs*rc_dep,1.05,color=BLUE,alpha=0.055,zorder=0)
a2.scatter(d.rg_vs_reference,d.rg_contaminant,s=110,color=BLUE,edgecolor="white",lw=1.5,zorder=3)
OFF2={"1 endorsed measure":(11,-3,"left"),"2 endorsed measures":(9,4,"left"),
      "3 endorsed measures":(9,-2,"left"),"4-5 endorsed measures":(9,4,"left"),
      "CIDI lifetime depression\n(reference)":(-8,-24,"right")}
for _,r in d.iterrows():
    dx,dy,ha=OFF2.get(r.label,(0,12,"center"))
    a2.annotate(r.label,(r.rg_vs_reference,r.rg_contaminant),textcoords="offset points",
                xytext=(dx,dy),ha=ha,fontsize=8,color=INK,linespacing=1.25)
a2.set_title("B. Depression\ncontamination is a property of the instrument",fontsize=10.6,loc="left",pad=10,weight="bold")
a2.set_xlabel("Genetic correlation with the reference definition\n(CIDI lifetime depression)",fontsize=8.6,color=MUTED)
a2.set_ylabel("Genetic correlation with the contaminant\n(neuroticism)",fontsize=8.6,color=MUTED)
a2.set_xlim(0.44,1.12); a2.set_ylim(0.42,0.95)

for ax,txt in ((a1,"above the line = contaminated"),(a2,"above the line = contaminated")):
    ax.grid(color=GRID,lw=0.7); ax.set_axisbelow(True)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
    ax.text(0.03,0.965,txt,transform=ax.transAxes,fontsize=8,color=MUTED,style="italic")
fig.suptitle("The same index, two trait families, two different failure modes",
             fontsize=13,color=INK,x=0.005,ha="left",y=1.075,weight="bold")
fig.text(0.005,1.005,"The dashed line is where a definition would sit if it were nothing but a noisier measure of the reference. "
         "Distance above it is contamination.",fontsize=8.3,color=MUTED,ha="left")
fig.savefig(f"{O}/Fig5_crossfamily.png",bbox_inches="tight",facecolor="white")
fig.savefig(f"{O}/Fig5_crossfamily.pdf",bbox_inches="tight",facecolor="white")
print("written")
