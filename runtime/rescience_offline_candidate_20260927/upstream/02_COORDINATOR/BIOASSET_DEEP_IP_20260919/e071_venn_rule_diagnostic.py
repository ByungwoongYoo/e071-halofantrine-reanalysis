#!/usr/bin/env python3
import io,gzip,re,json,math
import pandas as pd, numpy as np, requests

SERIES={"CEP290_LCA":("GSE327431","human"),"P23H":("GSE327432","rat"),"rd10":("GSE327434","mouse")}
FC=math.log2(1.2)
SOURCE={"H_only":470,"R_only":379,"M_only":405,"H_R_only":16,"R_M_only":16,"H_M_only":13,"all3":1}
ITPR3={"human":"ENSG00000096433","rat":"ENSRNOG00000052795","mouse":"ENSMUSG00000042644"}

def base(g):return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{g[:-3]}nnn/{g}/suppl/"
def get_deg(g):
 t=requests.get(base(g),timeout=120).text
 fs=[x for x in re.findall(r'href="([^"]+)"',t,re.I) if x.endswith(".csv.gz") and "deg" in x.lower()]
 r=requests.get(base(g)+fs[0],timeout=180);r.raise_for_status()
 d=pd.read_csv(gzip.GzipFile(fileobj=io.BytesIO(r.content)))
 d=d.copy();d["gene"]=d.iloc[:,0].astype(str).str.replace(r"\.\d+$","",regex=True)
 for c in ["logFC","P.Value","adj.P.Val","t"]:d[c]=pd.to_numeric(d[c],errors="coerce")
 return d.dropna(subset=["gene","logFC","P.Value","t"]).drop_duplicates("gene")

def biomart_all():
 xml='''<?xml version="1.0" encoding="UTF-8"?>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
<Dataset name="hsapiens_gene_ensembl" interface="default">
<Attribute name="ensembl_gene_id"/>
<Attribute name="mmusculus_homolog_ensembl_gene"/>
<Attribute name="rnorvegicus_homolog_ensembl_gene"/>
</Dataset></Query>'''
 for ep in ["https://www.ensembl.org/biomart/martservice","https://useast.ensembl.org/biomart/martservice"]:
  try:
   r=requests.get(ep,params={"query":xml},timeout=180);r.raise_for_status()
   z=pd.read_csv(io.StringIO(r.text),sep="\t")
   if len(z)>1000:break
  except Exception:z=None
 if z is None:raise RuntimeError("biomart")
 z.columns=["human","mouse","rat"]
 for c in z.columns:z[c]=z[c].fillna("").astype(str).str.replace(r"\.\d+$","",regex=True).str.strip()
 return z

D={m:get_deg(g) for m,(g,s) in SERIES.items()}
BM=biomart_all()
map_rat=BM[BM.rat!=""][["rat","human"]].drop_duplicates()
map_mouse=BM[BM.mouse!=""][["mouse","human"]].drop_duplicates()

def select(d,scheme):
 x=d.copy()
 if scheme=="top500_absfc":return x.nlargest(500,"logFC",keep="all") if False else x.assign(A=x.logFC.abs()).nlargest(500,"A")
 if scheme=="top500_p":return x.nsmallest(500,"P.Value")
 if scheme=="top500_abst":return x.assign(A=x.t.abs()).nlargest(500,"A")
 if scheme=="fc_top500_absfc":return x[x.logFC.abs()>=FC].assign(A=lambda q:q.logFC.abs()).nlargest(500,"A")
 if scheme=="p_top500_absfc":return x[x["P.Value"]<.05].assign(A=lambda q:q.logFC.abs()).nlargest(500,"A")
 if scheme=="pfc_top500_absfc":return x[(x["P.Value"]<.05)&(x.logFC.abs()>=FC)].assign(A=lambda q:q.logFC.abs()).nlargest(500,"A")
 if scheme=="top250_each":
  return pd.concat([x[x.logFC>0].nlargest(250,"logFC"),x[x.logFC<0].nsmallest(250,"logFC")]).drop_duplicates("gene")
 if scheme=="p_top250_each":
  q=x[x["P.Value"]<.05]; return pd.concat([q[q.logFC>0].nlargest(250,"logFC"),q[q.logFC<0].nsmallest(250,"logFC")]).drop_duplicates("gene")
 if scheme=="pfc_top250_each":
  q=x[(x["P.Value"]<.05)&(x.logFC.abs()>=FC)];return pd.concat([q[q.logFC>0].nlargest(250,"logFC"),q[q.logFC<0].nsmallest(250,"logFC")]).drop_duplicates("gene")
 raise KeyError(scheme)

def mapped_sets(scheme):
 h=set(select(D["CEP290_LCA"],scheme).gene)
 rsel=select(D["P23H"],scheme)[["gene"]].rename(columns={"gene":"rat"})
 msel=select(D["rd10"],scheme)[["gene"]].rename(columns={"gene":"mouse"})
 r=set(rsel.merge(map_rat,on="rat").human)
 m=set(msel.merge(map_mouse,on="mouse").human)
 return h,r,m

def regions(h,r,m):
 return {
  "H_only":len(h-r-m),
  "R_only":len(r-h-m),
  "M_only":len(m-h-r),
  "H_R_only":len((h&r)-m),
  "R_M_only":len((r&m)-h),
  "H_M_only":len((h&m)-r),
  "all3":len(h&r&m),
 }

schemes=["top500_absfc","top500_p","top500_abst","fc_top500_absfc","p_top500_absfc","pfc_top500_absfc","top250_each","p_top250_each","pfc_top250_each"]
out={}
for s in schemes:
 try:
  h,r,m=mapped_sets(s); reg=regions(h,r,m)
  score=sum(abs(reg[k]-SOURCE[k]) for k in SOURCE)
  out[s]={"pre_map":{"human":len(select(D["CEP290_LCA"],s)),"rat":len(select(D["P23H"],s)),"mouse":len(select(D["rd10"],s))},
          "mapped_sizes":{"human":len(h),"rat":len(r),"mouse":len(m)},"regions":reg,"L1_source_distance":score,
          "ITPR3_in":{"human":ITPR3["human"] in h,"rat":ITPR3["human"] in r,"mouse":ITPR3["human"] in m},
          "all3_ids":sorted(h&r&m)[:30]}
 except Exception as e: out[s]={"error":repr(e)}
print("E071_VENN_RULE_DIAGNOSTIC_JSON_BEGIN")
print(json.dumps({"source":SOURCE,"schemes":out},indent=2,sort_keys=True))
print("E071_VENN_RULE_DIAGNOSTIC_JSON_END")
