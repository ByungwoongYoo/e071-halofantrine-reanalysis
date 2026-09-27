#!/usr/bin/env python3
import io, gzip, json, math, re
import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr, pearsonr

SERIES = {
    "CEP290_LCA": ("GSE327431","human"),
    "P23H": ("GSE327432","rat"),
    "rd10": ("GSE327434","mouse"),
}
FC_THRESH = math.log2(1.2)
RNG = np.random.default_rng(20260921)

def directory_url(gse):
    prefix=gse[:-3]+"nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/suppl/"

def discover_deg(gse):
    u=directory_url(gse)
    r=requests.get(u,timeout=120); r.raise_for_status()
    names=re.findall(r'href="([^"]+)"',r.text,flags=re.I)
    cands=[x for x in names if x.lower().endswith(".csv.gz") and "deg" in x.lower()]
    if not cands:
        raise RuntimeError(f"{gse}: no DEG csv.gz in {u}; names={names[:100]}")
    print("DISCOVER",gse,cands)
    return u+cands[0]

def fetch_csv_gz(url):
    r=requests.get(url,timeout=180); r.raise_for_status()
    return pd.read_csv(gzip.GzipFile(fileobj=io.BytesIO(r.content)))

def col_exact(cols,names):
    lut={str(c).lower():c for c in cols}
    for n in names:
        if n.lower() in lut: return lut[n.lower()]
    return None

def standardize(df,model):
    cols=list(df.columns)
    gene=cols[0]
    lfc=col_exact(cols,["logFC","log2FoldChange","log2FC"])
    p=col_exact(cols,["P.Value","pvalue","p_value","PValue","P.Val"])
    adj=col_exact(cols,["adj.P.Val","FDR","padj","qvalue"])
    if lfc is None or p is None:
        raise RuntimeError(f"{model}: columns={cols}")
    z=pd.DataFrame({
        "gene_id":df[gene].astype(str).str.replace(r"\.\d+$","",regex=True).str.strip(),
        "logFC":pd.to_numeric(df[lfc],errors="coerce"),
        "p":pd.to_numeric(df[p],errors="coerce"),
    })
    if adj is not None:
        z["adj"]=pd.to_numeric(df[adj],errors="coerce")
    else:
        pv=np.nan_to_num(z["p"].to_numpy(float),nan=1.0)
        order=np.argsort(pv); ranked=pv[order]
        bh=ranked*len(ranked)/(np.arange(len(ranked))+1)
        bh=np.minimum.accumulate(bh[::-1])[::-1]
        q=np.empty(len(bh)); q[order]=np.minimum(bh,1.0); z["adj"]=q
    z=z.dropna(subset=["gene_id","logFC","p"]).copy()
    z=z[(z.gene_id!="")&(z.gene_id!="nan")]
    z=z.sort_values("p").drop_duplicates("gene_id")
    return z, {"gene_col":str(gene),"lfc_col":str(lfc),"p_col":str(p),"adj_col":str(adj)}

def author_list(z):
    q=z[(z.p<0.05)&(z.logFC.abs()>=FC_THRESH)].copy()
    up=q[q.logFC>0].nlargest(500,"logFC")
    dn=q[q.logFC<0].nsmallest(500,"logFC")
    return pd.concat([up,dn]).drop_duplicates("gene_id")

def fdr_list(z):
    return z[(z.adj<0.05)&(z.logFC.abs()>=FC_THRESH)].copy()

def pair_metrics(a,b,key="human"):
    m=a[[key,"logFC"]].merge(b[[key,"logFC"]],on=key,suffixes=("_a","_b"))
    x=m.logFC_a.to_numpy(float); y=m.logFC_b.to_numpy(float)
    ok=np.isfinite(x)&np.isfinite(y); x=x[ok]; y=y[ok]
    if len(x)<3: return {"n":int(len(x))}
    sign=float(np.mean(np.sign(x)==np.sign(y)))
    B=5000
    null=np.array([np.mean(np.sign(x)==np.sign(RNG.permutation(y))) for _ in range(B)])
    return {
      "n":int(len(x)),
      "spearman":float(spearmanr(x,y).statistic),
      "pearson":float(pearsonr(x,y).statistic),
      "cosine":float(np.dot(x,y)/(np.linalg.norm(x)*np.linalg.norm(y))),
      "sign_concordance":sign,
      "null_sign_mean":float(null.mean()),
      "sign_perm_p":float((1+(null>=sign).sum())/(B+1)),
    }

def overlap_stats(sets,universe):
    names=list(sets); arr=np.array(sorted(universe),dtype=object); U=len(arr)
    sizes=[len(sets[n]&universe) for n in names]
    obs2=sum(sum(g in sets[n] for n in names)>=2 for g in universe)
    obs3=len(set.intersection(*(sets[n]&universe for n in names)))
    B=10000
    n2=np.empty(B,dtype=int); n3=np.empty(B,dtype=int)
    for i in range(B):
        s=[set(RNG.choice(arr,size=k,replace=False)) for k in sizes]
        uu=set().union(*s)
        n2[i]=sum(sum(g in x for x in s)>=2 for g in uu)
        n3[i]=len(s[0]&s[1]&s[2])
    return {
      "universe":U,
      "sizes":dict(zip(names,map(int,sizes))),
      "observed_ge2":int(obs2),"observed_all3":int(obs3),
      "null_ge2_mean":float(n2.mean()),"null_ge2_sd":float(n2.std(ddof=1)),
      "null_all3_mean":float(n3.mean()),"null_all3_sd":float(n3.std(ddof=1)),
      "p_ge2":float((1+(n2>=obs2).sum())/(B+1)),
      "p_all3":float((1+(n3>=obs3).sum())/(B+1)),
    }

def biomart():
    from pathlib import Path
    p = Path("02_COORDINATOR/BIOASSET_DEEP_IP_20260919/SREP_SUBMISSION_PACKAGE_V1/ortholog_mapping/ensembl_one2one_core.tsv")
    if not p.is_file():
        raise FileNotFoundError(f"Frozen core mapping missing: {p}")
    t = pd.read_csv(p, sep="\\t", dtype=str).fillna("")
    required = ["human", "mouse", "rat"]
    if list(t.columns) != required:
        raise RuntimeError(f"Unexpected frozen core mapping columns: {list(t.columns)}")
    if len(t) != 15956:
        raise RuntimeError(f"Unexpected frozen core mapping row count: {len(t)}")
    if any(t[c].duplicated().any() for c in required):
        raise RuntimeError("Frozen core mapping violates unique one-to-one contract")
    print("BIOMART_FROZEN_CORE", len(t), str(p))
    return t[required].copy()


print("E071_GSE_REPLAY_START")
std={}; author={}; fdr={}; schema={}; urls={}
for name,(gse,species) in SERIES.items():
    u=discover_deg(gse); urls[name]=u
    print("DOWNLOAD",name,u)
    d=fetch_csv_gz(u)
    print("RAW",name,"shape",d.shape,"columns",list(map(str,d.columns)))
    print("RAW_HEAD",name,json.dumps(d.head(3).to_dict("records"),default=str))
    z,m=standardize(d,name); std[name]=z; schema[name]=m
    author[name]=author_list(z); fdr[name]=fdr_list(z)
    print("MODEL",name,"tested",len(z),"nominal_fc",len(z[(z.p<.05)&(z.logFC.abs()>=FC_THRESH)]),
          "author",len(author[name]),"fdr",len(fdr[name]),
          "fdr_min",float(z.adj.min()))

res={"urls":urls,"schema":schema,"models":{}}
for n in SERIES:
    z=std[n]; a=author[n]; q=fdr[n]
    res["models"][n]={
      "tested":int(len(z)),
      "nominal_p05_fc12":int(len(z[(z.p<.05)&(z.logFC.abs()>=FC_THRESH)])),
      "author_selected":int(len(a)),
      "author_up":int((a.logFC>0).sum()),"author_down":int((a.logFC<0).sum()),
      "fdr_selected":int(len(q)),
      "fdr_up":int((q.logFC>0).sum()),"fdr_down":int((q.logFC<0).sum()),
      "min_fdr":float(z.adj.min()),
    }

try:
    bm=biomart(); print("BIOMART_TRIPLES",len(bm))
    key={"CEP290_LCA":"human","rd10":"mouse","P23H":"rat"}
    allf={}; af={}; qf={}
    for n in SERIES:
        k=key[n]
        allf[n]=std[n].merge(bm,left_on="gene_id",right_on=k).drop_duplicates("human")
        af[n]=author[n].merge(bm,left_on="gene_id",right_on=k).drop_duplicates("human")
        qf[n]=fdr[n].merge(bm,left_on="gene_id",right_on=k).drop_duplicates("human")
    universe=set.intersection(*(set(allf[n].human) for n in SERIES))
    res["one2one_common_universe"]=len(universe)
    res["all_gene_pair_metrics"]={}
    for x,y in [("CEP290_LCA","P23H"),("CEP290_LCA","rd10"),("P23H","rd10")]:
        res["all_gene_pair_metrics"][f"{x}__{y}"]=pair_metrics(allf[x],allf[y])
    aset={n:set(af[n].human) for n in SERIES}
    qset={n:set(qf[n].human) for n in SERIES}
    res["author_overlap"]=overlap_stats(aset,universe)
    res["fdr_overlap"]=overlap_stats(qset,universe)
    common3=set.intersection(*(aset[n] for n in SERIES))
    common2=set(g for g in universe if sum(g in aset[n] for n in SERIES)>=2)
    res["author_common_all3_ids"]=sorted(common3)
    res["author_common_ge2_count"]=len(common2)
    # direction consistency for genes selected in >=2 models
    direction={}
    for h in common2:
        vals={}
        for n in SERIES:
            sub=af[n][af[n].human==h]
            if len(sub): vals[n]=float(sub.logFC.iloc[0])
        direction[h]=vals
    res["author_common_ge2_direction_consistent"]=sum(
        len(v)>=2 and len(set(np.sign(list(v.values()))))==1 for v in direction.values())
    res["author_common_ge2_direction_inconsistent"]=sum(
        len(v)>=2 and len(set(np.sign(list(v.values()))))>1 for v in direction.values())
except Exception as e:
    res["biomart_error"]=repr(e); print("BIOMART_ERROR",repr(e))

print("E071_RESULT_JSON_BEGIN")
print(json.dumps(res,indent=2,sort_keys=True))
print("E071_RESULT_JSON_END")
