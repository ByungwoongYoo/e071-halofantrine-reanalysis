#!/usr/bin/env python3
import io,gzip,json,math,re,itertools
import numpy as np, pandas as pd, requests
from scipy.stats import spearmanr, pearsonr

SERIES={
 "CEP290_LCA":("GSE327431","human"),
 "P23H":("GSE327432","rat"),
 "rd10":("GSE327434","mouse")
}
PSEUDO=0.5
FC=math.log2(1.2)

def durl(gse): return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/"
def files(gse):
 t=requests.get(durl(gse),timeout=120).text
 return [x for x in re.findall(r'href="([^"]+)"',t,re.I) if x.endswith(".gz")]
def fetch(gse,needle):
 fs=[x for x in files(gse) if all(k in x.lower() for k in needle)]
 if not fs: raise RuntimeError((gse,needle,files(gse)))
 r=requests.get(durl(gse)+fs[0],timeout=180); r.raise_for_status()
 return fs[0],pd.read_csv(gzip.GzipFile(fileobj=io.BytesIO(r.content)))
def groups(model,cols):
 if model=="CEP290_LCA":
  return [c for c in cols if c.startswith("PEN1B.HAL.")],[c for c in cols if c.startswith("PEN1B.DMSO.")]
 if model=="P23H":
  return [c for c in cols if c.startswith("P23H.F.HAL.")],[c for c in cols if c.startswith("P23H.F.DMSO.")]
 if model=="rd10":
  return [c for c in cols if c.startswith("rd10.F.HAL.")],[c for c in cols if c.startswith("rd10.F.DMSO.")]
 raise KeyError(model)
def vec(logx,H,D): return logx[H].mean(axis=1)-logx[D].mean(axis=1)
def corr(x,y):
 ok=np.isfinite(x)&np.isfinite(y)
 x=x[ok];y=y[ok]
 return {"n":int(len(x)),"pearson":float(pearsonr(x,y).statistic),"spearman":float(spearmanr(x,y).statistic),
         "sign":float(np.mean(np.sign(x)==np.sign(y)))}
def exact_centroid_perm(logx,H,D):
 cols=H+D;nH=len(H)
 X=logx[cols].to_numpy(float)
 sd=np.nanstd(X,axis=1,ddof=1); ok=np.isfinite(sd)&(sd>0)
 X=X[ok]; sd=sd[ok]
 # z-score per gene across disease samples; removes expression-scale dominance
 X=(X-np.nanmean(X,axis=1,keepdims=True))/sd[:,None]
 obs_idx=[cols.index(c) for c in H]
 def score(idx):
  idx=set(idx); a=[i for i in range(len(cols)) if i in idx]; b=[i for i in range(len(cols)) if i not in idx]
  d=np.nanmean(X[:,a],axis=1)-np.nanmean(X[:,b],axis=1)
  return float(np.nanmean(d*d))
 obs=score(obs_idx)
 vals=np.array([score(comb) for comb in itertools.combinations(range(len(cols)),nH)])
 return {"n_permutations":int(len(vals)),"score":obs,"null_mean":float(vals.mean()),"null_sd":float(vals.std(ddof=1)),
         "exact_p_ge":float(np.mean(vals>=obs)),"percentile":float(np.mean(vals<=obs))}

def biomart():
 xml='''<?xml version="1.0" encoding="UTF-8"?>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="hsapiens_gene_ensembl" interface="default">
<Attribute name="ensembl_gene_id"/>
<Attribute name="mmusculus_homolog_ensembl_gene"/>
<Attribute name="mmusculus_homolog_orthology_type"/>
<Attribute name="rnorvegicus_homolog_ensembl_gene"/>
<Attribute name="rnorvegicus_homolog_orthology_type"/>
</Dataset></Query>'''
 z=None
 for ep in ["https://www.ensembl.org/biomart/martservice","https://useast.ensembl.org/biomart/martservice"]:
  try:
   r=requests.get(ep,params={"query":xml},timeout=180);r.raise_for_status()
   q=pd.read_csv(io.StringIO(r.text),sep="\t")
   if len(q)>1000: z=q;break
  except Exception: pass
 if z is None: raise RuntimeError("biomart")
 z.columns=["human","mouse","mouse_type","rat","rat_type"]
 z=z[(z.mouse_type=="ortholog_one2one")&(z.rat_type=="ortholog_one2one")].copy()
 for c in ["human","mouse","rat"]:
  z[c]=z[c].fillna("").astype(str).str.replace(r"\.\d+$","",regex=True).str.strip()
 z=z[(z.human!="")&(z.mouse!="")&(z.rat!="")]
 for c in ["human","mouse","rat"]:
  vc=z[c].value_counts(); z=z[z[c].isin(vc[vc==1].index)]
 return z[["human","mouse","rat"]].drop_duplicates()

models={}; out={"models":{}}
for model,(gse,species) in SERIES.items():
 fn,cpm=fetch(gse,["gene","cpm","normalized"])
 dfn,deg=fetch(gse,["deg"])
 gene=cpm.iloc[:,0].astype(str).str.replace(r"\.\d+$","",regex=True)
 cpm=cpm.copy();cpm["_gene"]=gene
 sample=[c for c in cpm.columns if c not in [cpm.columns[0],"_gene"]]
 H,D=groups(model,sample)
 logx=np.log2(cpm[sample].astype(float)+PSEUDO)
 logx.index=gene
 full=vec(logx,H,D)
 full.name="sample_logFC"
 # validate against deposited logFC
 dg=pd.DataFrame({"gene":deg.iloc[:,0].astype(str).str.replace(r"\.\d+$","",regex=True),
                  "dep_logFC":pd.to_numeric(deg["logFC"],errors="coerce")}).dropna().drop_duplicates("gene")
 cmp=pd.DataFrame({"gene":full.index,"sample_logFC":full.values}).merge(dg,on="gene")
 validation=corr(cmp.sample_logFC.to_numpy(),cmp.dep_logFC.to_numpy())
 # LOO one sample at a time
 loo=[]
 for drop in H+D:
  h=[x for x in H if x!=drop]; d=[x for x in D if x!=drop]
  v=vec(logx,h,d).reindex(full.index)
  allc=corr(full.to_numpy(),v.to_numpy())
  hi=np.abs(full.to_numpy())>=FC
  hic=corr(full.to_numpy()[hi],v.to_numpy()[hi]) if hi.sum()>2 else {}
  loo.append({"drop":drop,"group":"HAL" if drop in H else "DMSO","all":allc,"high_effect":hic})
 # exact disease-sample label permutation separation
 perm=exact_centroid_perm(logx,H,D)
 out["models"][model]={
  "cpm_file":fn,"deg_file":dfn,"HAL":H,"DMSO":D,
  "validation_sample_vs_deposited":validation,
  "exact_label_permutation":perm,
  "loo":loo,
  "loo_all_spearman_min":float(min(x["all"]["spearman"] for x in loo)),
  "loo_all_spearman_median":float(np.median([x["all"]["spearman"] for x in loo])),
  "loo_high_effect_spearman_min":float(min(x["high_effect"]["spearman"] for x in loo)),
  "loo_high_effect_spearman_median":float(np.median([x["high_effect"]["spearman"] for x in loo])),
 }
 models[model]={"gene":full.index.to_series().astype(str).str.replace(r"\.\d+$","",regex=True).values,
                "full":full.to_numpy(),"logx":logx,"H":H,"D":D,"loo":loo}

# construct full + LOO vectors by model
states={}
for model in SERIES:
 logx=models[model]["logx"]; H=models[model]["H"];D=models[model]["D"]
 ss={"FULL":vec(logx,H,D)}
 for drop in H+D:
  ss["DROP:"+drop]=vec(logx,[x for x in H if x!=drop],[x for x in D if x!=drop])
 states[model]=ss

bm=biomart(); key={"CEP290_LCA":"human","P23H":"rat","rd10":"mouse"}
mapped={}
for model in SERIES:
 base=pd.DataFrame({"gene":states[model]["FULL"].index})
 base["gene"]=base.gene.str.replace(r"\.\d+$","",regex=True)
 mp=base.merge(bm,left_on="gene",right_on=key[model])[["gene","human"]].drop_duplicates("human")
 mapped[model]=mp
out["cross_model_loo"]={}
for a,b in [("CEP290_LCA","P23H"),("CEP290_LCA","rd10"),("P23H","rd10")]:
 vals=[]
 for sa,va in states[a].items():
  da=pd.DataFrame({"gene":va.index.str.replace(r"\.\d+$","",regex=True),"x":va.values}).merge(mapped[a],on="gene").drop_duplicates("human")
  for sb,vb in states[b].items():
   db=pd.DataFrame({"gene":vb.index.str.replace(r"\.\d+$","",regex=True),"y":vb.values}).merge(mapped[b],on="gene").drop_duplicates("human")
   m=da[["human","x"]].merge(db[["human","y"]],on="human").dropna()
   rho=float(spearmanr(m.x,m.y).statistic)
   vals.append({"a":sa,"b":sb,"n":int(len(m)),"spearman":rho})
 arr=np.array([x["spearman"] for x in vals])
 full=[x for x in vals if x["a"]=="FULL" and x["b"]=="FULL"][0]
 out["cross_model_loo"][a+"__"+b]={
  "full":full,"all_state_pairs":len(vals),"min":float(arr.min()),"median":float(np.median(arr)),"max":float(arr.max()),
  "negative_fraction":float(np.mean(arr<0))
 }

print("E071_SAMPLE_LEVEL_ROBUSTNESS_JSON_BEGIN")
print(json.dumps(out,indent=2,sort_keys=True))
print("E071_SAMPLE_LEVEL_ROBUSTNESS_JSON_END")
