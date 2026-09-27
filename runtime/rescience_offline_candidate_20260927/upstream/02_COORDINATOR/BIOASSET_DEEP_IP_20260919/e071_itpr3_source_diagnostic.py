#!/usr/bin/env python3
import io,gzip,re,json,math
import numpy as np,pandas as pd,requests
from scipy.stats import ttest_ind, mannwhitneyu

CFG={
 "CEP290_LCA":{"gse":"GSE327431","id":"ENSG00000096433","hal":"PEN1B.HAL.","dmso":"PEN1B.DMSO."},
 "P23H":{"gse":"GSE327432","id":"ENSRNOG00000052795","hal":"P23H.F.HAL.","dmso":"P23H.F.DMSO."},
 "rd10":{"gse":"GSE327434","id":"ENSMUSG00000042644","hal":"rd10.F.HAL.","dmso":"rd10.F.DMSO."},
}
def base(g):return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/suppl/"
def ls(g):
 t=requests.get(base(g),timeout=120).text
 return [x for x in re.findall(r'href="([^"]+)"',t,re.I) if x.endswith(".gz")]
def get(g,preds):
 fs=[x for x in ls(g) if all(q in x.lower() for q in preds)]
 if not fs:raise RuntimeError((g,preds))
 r=requests.get(base(g)+fs[0],timeout=180);r.raise_for_status()
 return fs[0],pd.read_csv(gzip.GzipFile(fileobj=io.BytesIO(r.content)))
out={}
for model,cfg in CFG.items():
 dfn,d=get(cfg["gse"],["deg"])
 gene=d.iloc[:,0].astype(str).str.replace(r"\.\d+$","",regex=True)
 d=d.copy();d["_gene"]=gene
 row=d[d._gene==cfg["id"]].iloc[0]
 p=pd.to_numeric(d["P.Value"],errors="coerce")
 lfc=pd.to_numeric(d["logFC"],errors="coerce")
 t=pd.to_numeric(d["t"],errors="coerce")
 # 1-based ranks
 rank_p=int(p.rank(method="min",ascending=True).loc[row.name])
 rank_abs_t=int(t.abs().rank(method="min",ascending=False).loc[row.name])
 rank_abs_fc=int(lfc.abs().rank(method="min",ascending=False).loc[row.name])
 positive=d[lfc>0].copy(); positive["_lfc"]=pd.to_numeric(positive["logFC"],errors="coerce")
 negative=d[lfc<0].copy(); negative["_lfc"]=pd.to_numeric(negative["logFC"],errors="coerce")
 if float(row["logFC"])>0:
  direction_rank=int(positive["_lfc"].rank(method="min",ascending=False).loc[row.name])
 else:
  direction_rank=int((-negative["_lfc"]).rank(method="min",ascending=False).loc[row.name])

 cfn,cpm=get(cfg["gse"],["gene","cpm","normalized"])
 cpm=cpm.copy(); cpm["_gene"]=cpm.iloc[:,0].astype(str).str.replace(r"\.\d+$","",regex=True)
 cr=cpm[cpm._gene==cfg["id"]].iloc[0]
 H=[c for c in cpm.columns if str(c).startswith(cfg["hal"])]
 D=[c for c in cpm.columns if str(c).startswith(cfg["dmso"])]
 hv=np.log2(pd.to_numeric(cr[H],errors="coerce").to_numpy(float)+0.5)
 dv=np.log2(pd.to_numeric(cr[D],errors="coerce").to_numpy(float)+0.5)
 welch=ttest_ind(hv,dv,equal_var=False)
 mw=mannwhitneyu(hv,dv,alternative="two-sided",method="exact")
 out[model]={
  "deg_file":dfn,"cpm_file":cfn,
  "logFC":float(row["logFC"]),"p":float(row["P.Value"]),"fdr":float(row["adj.P.Val"]),"t":float(row["t"]),
  "rank_p":rank_p,"rank_abs_t":rank_abs_t,"rank_abs_logFC":rank_abs_fc,"direction_rank_by_logFC":direction_rank,
  "passes_p05":bool(float(row["P.Value"])<0.05),
  "passes_fc12":bool(abs(float(row["logFC"]))>=math.log2(1.2)),
  "HAL_n":len(H),"DMSO_n":len(D),
  "HAL_log2cpm_mean":float(hv.mean()),"DMSO_log2cpm_mean":float(dv.mean()),
  "sample_log2_difference":float(hv.mean()-dv.mean()),
  "welch_t_p":float(welch.pvalue),"mannwhitney_exact_p":float(mw.pvalue),
  "HAL_samples":H,"DMSO_samples":D
 }
print("E071_ITPR3_SOURCE_DIAGNOSTIC_JSON_BEGIN")
print(json.dumps(out,indent=2,sort_keys=True))
print("E071_ITPR3_SOURCE_DIAGNOSTIC_JSON_END")
