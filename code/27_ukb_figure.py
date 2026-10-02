import os
ROOT=os.environ.get("ROOT","/tmp")
import pandas as pd, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
T=pd.read_csv(f'{ROOT}/ukb/Table8_ukb_definition_scan.csv')
A=T[~T.duplicate_definition].sort_values('excess_max')
COL={'self-report':'#c44e52','touchscreen composite':'#dd8452','curated endpoint':'#4c72b0',
     'doctor-diagnosed Q':'#55a868','hospital ICD':'#8172b3'}
fig,axes=plt.subplots(1,2,figsize=(14,8.6),gridspec_kw={'width_ratios':[2.25,1]})
ax=axes[0]; y=np.arange(len(A))
ax.barh(y,A.excess_max,color=[COL.get(t,'#888') for t in A.def_type],height=0.72,zorder=3)
ax.errorbar(A.excess_max,y,xerr=A.excess_max_se,fmt='none',ecolor='0.25',elinewidth=0.8,capsize=1.5,zorder=4)
lab=[f"{f}   ({d} vs {r})" for f,d,r in zip(A.family,A.def_type,A.ref_type)]
ax.set_yticks(y); ax.set_yticklabels(lab,fontsize=7.5)
ax.axvline(0,color='k',lw=0.8,zorder=5)
ax.axvline(0.376,color='#c44e52',ls='--',lw=1.1,zorder=5)
ax.text(0.376,len(A)-0.2,' AD self-report (+0.38)',color='#c44e52',fontsize=7.5,va='top')
dup=T[T.duplicate_definition]
ax.set_xlabel('Excess genetic correlation with worst-case contaminant of the panel\n'
              r'excess = $r_g(D,C)-r_g(D,R)\,r_g(R,C)$',fontsize=9)
ax.set_title(f'Phenotype-definition contamination across {T.family.nunique()} UK Biobank trait families\n'
             f'one cohort, one array, one control set; {len(A)} definition-vs-reference pairs, '
             f'{len(dup)} a-priori negative controls (curated endpoint = same ICD code): max excess {dup.excess_max.abs().max():.3f}',fontsize=10,loc='left')
used=[t for t in COL if t in set(A.def_type)]
h=[plt.Rectangle((0,0),1,1,color=COL[t]) for t in used]
ax.legend(h,used,fontsize=7.5,loc='lower right',title='definition under test',title_fontsize=8)
ax.grid(axis='x',color='0.9',zorder=0)
ax=axes[1]
o=A.groupby('def_type').excess_max.median().sort_values()
rng=np.random.default_rng(0)
for i,(t,v) in enumerate(o.items()):
    s=A[A.def_type==t]
    ax.scatter(s.excess_max,np.full(len(s),i)+rng.normal(0,.07,len(s)),s=30,
               color=COL.get(t,'#888'),alpha=.85,edgecolor='w',linewidth=.5,zorder=3)
    ax.plot([v,v],[i-.28,i+.28],color='k',lw=2.2,zorder=4)
ax.scatter(dup.excess_max,np.full(len(dup),-1)+rng.normal(0,.07,len(dup)),s=30,color='0.6',
           edgecolor='w',linewidth=.5,zorder=3)
ax.set_yticks(list(range(-1,len(o)))); ax.set_yticklabels(['duplicate definition\n(negative control)']+list(o.index),fontsize=8)
ax.axvline(0,color='k',lw=0.8); ax.set_xlabel('Excess $r_g$',fontsize=9)
ax.set_title('By type of case definition (bar = median)',fontsize=10,loc='left',pad=26)
ax.grid(axis='x',color='0.9',zorder=0)
for a in axes: a.spines[['top','right']].set_visible(False)
plt.tight_layout()
plt.savefig(f'{ROOT}/ukb/Fig6_ukb_definition_scan.png',dpi=200,bbox_inches='tight')
plt.savefig(f'{ROOT}/ukb/Fig6_ukb_definition_scan.pdf',bbox_inches='tight')
print("ok",len(A),len(dup))
