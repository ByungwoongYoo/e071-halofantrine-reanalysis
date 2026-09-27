#!/usr/bin/env python3
import io, gzip, json, math, re, tempfile
import numpy as np
import pandas as pd
import requests
import gseapy as gp

SERIES = {
    "CEP290_LCA": ("GSE327431","human"),
    "P23H": ("GSE327432","rat"),
    "rd10": ("GSE327434","mouse"),
}

def directory_url(gse):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse[:-3]}nnn/{gse}/suppl/"

def discover_deg(gse):
    u=directory_url(gse)
    t=requests.get(u,timeout=120).text
    c=[x for x in re.findall(r'href="([^"]+)"',t,re.I) if x.lower().endswith(".csv.gz") and "deg" in x.lower()]
    if not c: raise RuntimeError((gse,u))
    return u+c[0]

def fetch(url):
    r=requests.get(url,timeout=180); r.raise_for_status()
    return pd.read_csv(gzip.GzipFile(fileobj=io.BytesIO(r.content)))

def biomart():
    from pathlib import Path
    p = Path("02_COORDINATOR/BIOASSET_DEEP_IP_20260919/SREP_SUBMISSION_PACKAGE_V1/ortholog_mapping/ensembl_one2one_gsea.tsv")
    if not p.is_file():
        raise FileNotFoundError(f"Frozen GSEA mapping missing: {p}")
    z = pd.read_csv(p, sep="\\t", dtype=str).fillna("")
    required = ["human", "symbol", "mouse", "rat"]
    if list(z.columns) != required:
        raise RuntimeError(f"Unexpected frozen GSEA mapping columns: {list(z.columns)}")
    if len(z) != 15940:
        raise RuntimeError(f"Unexpected frozen GSEA mapping row count: {len(z)}")
    for c in ["human", "mouse", "rat"]:
        if z[c].duplicated().any():
            raise RuntimeError(f"Frozen GSEA mapping violates unique one-to-one contract: {c}")
    print("BIOMART_FROZEN_GSEA", len(z), str(p))
    return z[required].copy()


bm=biomart()
print("BIOMART",len(bm))
key={"CEP290_LCA":"human","rd10":"mouse","P23H":"rat"}
rankings={}
for name,(gse,species) in SERIES.items():
    d=fetch(discover_deg(gse))
    gene=d.iloc[:,0].astype(str).str.replace(r"\.\d+$","",regex=True)
    x=pd.DataFrame({"gene_id":gene,"logFC":pd.to_numeric(d["logFC"],errors="coerce"),
                    "p":pd.to_numeric(d["P.Value"],errors="coerce"),
                    "t":pd.to_numeric(d["t"],errors="coerce")})
    x=x.dropna()
    x=x.merge(bm,left_on="gene_id",right_on=key[name]).drop_duplicates("symbol")
    # use source moderated t statistic as threshold-free ranking; deterministic tie-break by symbol.
    x=x.sort_values(["t","symbol"],ascending=[False,True])
    rankings[name]=x[["symbol","t","logFC","p"]].copy()
    print("RANK",name,len(x),"t range",float(x.t.min()),float(x.t.max()))

# fixed common symbol universe only, so all models test identical genes.
common=set.intersection(*(set(x.symbol) for x in rankings.values()))
print("COMMON_SYMBOL_UNIVERSE",len(common))
for n in rankings:
    rankings[n]=rankings[n][rankings[n].symbol.isin(common)].drop_duplicates("symbol")

lib=gp.get_library(name="KEGG_2021_Human",organism="Human")
# Restrict pathways to genes in common universe and stable sizes.
lib2={}
for term,genes in lib.items():
    gs=[g for g in genes if g in common]
    if 10 <= len(gs) <= 500:
        lib2[term]=gs
print("KEGG_PATHWAYS",len(lib2))

results={}
for name,r in rankings.items():
    print("GSEA_START",name)
    pre=gp.prerank(
        rnk=r[["symbol","t"]],
        gene_sets=lib2,
        min_size=10,max_size=500,
        permutation_num=2000,
        seed=20260921,
        threads=4,
        verbose=False,
    )
    df=pre.res2d.copy()
    # normalize column names across gseapy versions
    cols={str(c).lower():c for c in df.columns}
    termcol=cols.get("term","Term")
    nescol=cols.get("nes","NES")
    fdrcol=cols.get("fdr q-val","FDR q-val")
    nomcol=cols.get("nom p-val","NOM p-val")
    z=df[[termcol,nescol,fdrcol,nomcol]].copy()
    z.columns=["term","NES","FDR","NOM_P"]
    z["NES"]=pd.to_numeric(z.NES,errors="coerce")
    z["FDR"]=pd.to_numeric(z.FDR,errors="coerce")
    z["NOM_P"]=pd.to_numeric(z.NOM_P,errors="coerce")
    z=z.dropna(subset=["NES","FDR"])
    results[name]=z
    print("GSEA_DONE",name,"FDR05",int((z.FDR<.05).sum()),"FDR25",int((z.FDR<.25).sum()))
    print("TOP",name,z.sort_values("FDR").head(12).to_json(orient="records"))

# Cross-model pathway NES geometry across all common KEGG terms with results in every model.
wide=None
for n,z in results.items():
    zz=z[["term","NES","FDR"]].rename(columns={"NES":f"NES_{n}","FDR":f"FDR_{n}"})
    wide=zz if wide is None else wide.merge(zz,on="term")
out={"common_gene_universe":len(common),"kegg_pathways":len(lib2),"models":{}}
for n,z in results.items():
    out["models"][n]={
      "FDR05_count":int((z.FDR<.05).sum()),
      "FDR25_count":int((z.FDR<.25).sum()),
      "FDR05_terms":z[z.FDR<.05].sort_values("FDR")[["term","NES","FDR"]].to_dict("records"),
    }

# terms significant at FDR<.05 in >=2/all3; and direction-consistent.
sigsets={n:set(results[n].loc[results[n].FDR<.05,"term"]) for n in results}
allterms=set(wide.term)
out["FDR05_ge2_terms"]=sorted([t for t in allterms if sum(t in sigsets[n] for n in results)>=2])
out["FDR05_all3_terms"]=sorted(set.intersection(*(sigsets[n] for n in results)))
# directional all-model NES sign consistency and correlations
nescols=[f"NES_{n}" for n in results]
W=wide.dropna(subset=nescols).copy()
out["all_pathway_NES_sign_all3_concordant"]=int((np.sign(W[nescols]).nunique(axis=1)==1).sum())
out["all_pathway_count"]=int(len(W))
out["pairwise_NES_corr"]={}
from scipy.stats import spearmanr, pearsonr
pairs=[("CEP290_LCA","P23H"),("CEP290_LCA","rd10"),("P23H","rd10")]
for a,b in pairs:
    x=W[f"NES_{a}"]; y=W[f"NES_{b}"]
    out["pairwise_NES_corr"][f"{a}__{b}"]={
      "spearman":float(spearmanr(x,y).statistic),
      "pearson":float(pearsonr(x,y).statistic),
      "sign_concordance":float(np.mean(np.sign(x)==np.sign(y))),
    }
# specifically track source-highlighted pathways by fuzzy terms
queries=["OXIDATIVE PHOSPHORYLATION","RIBOSOME","CIRCADIAN"]
out["highlighted"]={}
for q in queries:
    hit=W[W.term.str.upper().str.contains(q,regex=False)]
    out["highlighted"][q]=hit.to_dict("records")

print("E071_GSEA_JSON_BEGIN")
print(json.dumps(out,indent=2,sort_keys=True))
print("E071_GSEA_JSON_END")
