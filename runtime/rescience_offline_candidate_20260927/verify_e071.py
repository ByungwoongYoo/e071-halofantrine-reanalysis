"""Bounded, local, read-only-source E071 claim verification (2026-09-27).

Downloads public deposited data, preserves original compressed bytes, and does
not refit or claim to reproduce the source limma model or exact S7A pipeline.
"""
from pathlib import Path
import argparse
import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import json
import math
import platform
import re
import sys
import numpy as np
import pandas as pd
import requests
import scipy
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "source_data"
OUT = ROOT / "results"
OFFLINE = False
BASE = ROOT / "upstream/02_COORDINATOR/BIOASSET_DEEP_IP_20260919"
MAP = BASE / "SREP_SUBMISSION_PACKAGE_V1/ortholog_mapping/ensembl_one2one_core.tsv"
CFG = {
    "human": dict(gse="GSE327431", gene="ENSG00000096433", hal="PEN1B.HAL.", dmso="PEN1B.DMSO.", expected_n=[5, 5]),
    "rat": dict(gse="GSE327432", gene="ENSRNOG00000052795", hal="P23H.F.HAL.", dmso="P23H.F.DMSO.", expected_n=[4, 5]),
    "mouse": dict(gse="GSE327434", gene="ENSMUSG00000042644", hal="rd10.F.HAL.", dmso="rd10.F.DMSO.", expected_n=[3, 3]),
}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

def acquire(spec):
    species, cfg = spec
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/{cfg['gse']}/suppl/"
    index = DATA / f"{cfg['gse']}_suppl_index.html"
    if OFFLINE:
        listing = index.read_text(encoding="utf-8")
    else:
        r = requests.get(base, timeout=(20, 90)); r.raise_for_status()
        index.write_bytes(r.content)
        listing = r.text
    names = sorted(set(re.findall(r'href="([^"/]+\.csv\.gz)"', listing, re.I)))
    result = {}
    for role, keys in [("deg", ["deg"]), ("cpm", ["gene", "cpm", "normalized"])]:
        matches = [n for n in names if all(k in n.lower() for k in keys)]
        assert len(matches) == 1, (species, role, matches)
        name = matches[0]; dest = DATA / name
        if not dest.exists():
            assert not OFFLINE, f"Missing frozen input {dest}"
            z = requests.get(base + name, timeout=(20, 180)); z.raise_for_status()
            gzip.decompress(z.content)  # validates CRC and truncation before saving
            dest.write_bytes(z.content)
        gzip.decompress(dest.read_bytes())
        result[role] = {"filename": name, "url": base + name, "bytes": dest.stat().st_size, "sha256": sha(dest)}
    print("ACQUIRED", species, flush=True)
    return species, result

def load_table(path):
    d = pd.read_csv(path)
    assert d.iloc[:, 0].notna().all()
    d = d.rename(columns={d.columns[0]: "gene"})
    d["gene"] = d.gene.astype(str).str.replace(r"\.\d+$", "", regex=True)
    assert not d.gene.duplicated().any(), path
    return d

def select(d, rule):
    d = d.copy()
    d["absfc"] = d.logFC.abs(); d["abst"] = d.t.abs()
    if rule.startswith("pfc_"):
        d = d[(d["P.Value"] < .05) & (d.absfc >= math.log2(1.2))]
    elif rule.startswith("p_"):
        d = d[d["P.Value"] < .05]
    elif rule.startswith("fc_"):
        d = d[d.absfc >= math.log2(1.2)]
    if "top250_each" in rule:
        return set(pd.concat([d[d.logFC > 0].nlargest(250, "logFC"), d[d.logFC < 0].nsmallest(250, "logFC")]).gene)
    if rule == "top500_p":
        return set(d.nsmallest(500, "P.Value").gene)
    return set(d.nlargest(500, "abst" if rule == "top500_abst" else "absfc").gene)

def venn(a, b, c):
    return {"H_only":len(a-b-c),"R_only":len(b-a-c),"M_only":len(c-a-b),"H_R_only":len(a&b-c),"R_M_only":len(b&c-a),"H_M_only":len(a&c-b),"all3":len(a&b&c)}

def main():
    DATA.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
    sources = dict(cf.ThreadPoolExecutor(max_workers=3).map(acquire, CFG.items()))
    receipt = ROOT / "results/source_manifest.json"
    if OFFLINE:
        original = json.loads(receipt.read_text(encoding="utf-8"))
        assert sources == original["sources"], "Frozen source hashes/metadata changed"
        if OUT / "source_manifest.json" != receipt:
            (OUT / "source_manifest.json").write_bytes(receipt.read_bytes())
    else:
        save_json(OUT / "source_manifest.json", {"retrieved_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"sources":sources})
    tables = {}; samples = {}; diagnostics = {}; itpr3 = []; sample_rows = []
    for species, cfg in CFG.items():
        d = load_table(DATA / sources[species]["deg"]["filename"])
        c = load_table(DATA / sources[species]["cpm"]["filename"])
        for col in ["logFC", "P.Value", "adj.P.Val", "t"]:
            assert pd.api.types.is_numeric_dtype(d[col]) and np.isfinite(d[col]).all(), (species,col)
        assert d["P.Value"].between(0,1).all() and d["adj.P.Val"].between(0,1).all()
        assert set(d.gene) == set(c.gene)
        assert np.isfinite(c.iloc[:,1:].to_numpy()).all() and (c.iloc[:,1:].to_numpy() >= 0).all()
        h = [x for x in c if x.startswith(cfg["hal"])]; v = [x for x in c if x.startswith(cfg["dmso"])]
        assert [len(h),len(v)] == cfg["expected_n"] and not set(h)&set(v)
        y = c.set_index("gene").loc[d.gene]
        effects = np.log2(y[h]+.5).mean(axis=1)-np.log2(y[v]+.5).mean(axis=1)
        bh = stats.false_discovery_control(d["P.Value"].to_numpy(), method="bh")
        row = d.set_index("gene").loc[cfg["gene"]]
        cr = y.loc[cfg["gene"]]
        hv = np.log2(cr[h].to_numpy(float)+.5); vv = np.log2(cr[v].to_numpy(float)+.5)
        nominal = (d["P.Value"] < .05) & (d.logFC.abs() >= math.log2(1.2))
        fdr = (d["adj.P.Val"] < .05) & (d.logFC.abs() >= math.log2(1.2))
        itpr3.append(dict(species=species,gse=cfg["gse"],gene=cfg["gene"],comparison="disease HAL minus disease DMSO",logFC=float(row.logFC),nominal_p=float(row["P.Value"]),adjusted_p=float(row["adj.P.Val"]),passes_nominal_p05=bool(row["P.Value"]<.05),passes_adjusted_p05=bool(row["adj.P.Val"]<.05),passes_abs_fc12=bool(abs(row.logFC)>=math.log2(1.2)),abs_logFC_rank=int(d.logFC.abs().rank(method="min",ascending=False).loc[d.gene==cfg["gene"]].iloc[0]),p_rank=int(d["P.Value"].rank(method="min").loc[d.gene==cfg["gene"]].iloc[0]),HAL_n=len(h),DMSO_n=len(v),sample_log2_CPM_difference=float(hv.mean()-vv.mean()),welch_p=float(stats.ttest_ind(hv,vv,equal_var=False).pvalue),mannwhitney_two_sided_exact_p=float(stats.mannwhitneyu(hv,vv,alternative="two-sided",method="exact").pvalue)))
        diagnostics[species] = dict(tested_genes=len(d),nominal_p05_only=int((d["P.Value"]<.05).sum()),nominal_p05_fc12=int(nominal.sum()),fdr05_fc12=int(fdr.sum()),minimum_adjusted_p=float(d["adj.P.Val"].min()),BH_from_full_deposited_table_max_abs_difference=float(np.max(np.abs(bh-d["adj.P.Val"]))),CPM_vs_deposit_Pearson=float(stats.pearsonr(effects,d.logFC).statistic),CPM_vs_deposit_Spearman=float(stats.spearmanr(effects,d.logFC).statistic),CPM_vs_deposit_sign_agreement=float(np.mean(np.sign(effects.to_numpy())==np.sign(d.logFC.to_numpy()))),HAL_samples=h,DMSO_samples=v,excluded_CPM_samples=[x for x in c.columns[1:] if x not in h+v])
        for group, names in [("HAL",h),("DMSO",v)]:
            for name in names:
                sample_rows.append(dict(species=species,gene=cfg["gene"],group=group,sample=name,CPM=float(cr[name]),log2_CPM_plus_half=float(np.log2(cr[name]+.5))))
        tables[species]=d; samples[species]=effects
    pd.DataFrame(itpr3).to_csv(OUT / "itpr3_claim_verification.csv",index=False)
    pd.DataFrame(sample_rows).to_csv(OUT / "itpr3_sample_values.csv",index=False)
    orth = pd.read_csv(MAP,sep="\t",dtype=str)
    assert len(orth)==15956 and list(orth.columns)==["human","mouse","rat"]
    assert all(not orth[s].duplicated().any() for s in CFG)
    orth_itpr3 = orth[orth.human==CFG['human']['gene']]
    assert len(orth_itpr3)==1 and all(orth_itpr3.iloc[0][s]==CFG[s]['gene'] for s in CFG)
    pairwise=[]
    for a,b in [("human","rat"),("human","mouse"),("rat","mouse")]:
        p=orth[[a,b]].merge(tables[a][["gene","logFC"]],left_on=a,right_on="gene").drop(columns="gene").rename(columns={"logFC":"a"}).merge(tables[b][["gene","logFC"]],left_on=b,right_on="gene").rename(columns={"logFC":"b"})
        pairwise.append(dict(pair=a+"_"+b,n=len(p),spearman=float(stats.spearmanr(p.a,p.b).statistic),pearson=float(stats.pearsonr(p.a,p.b).statistic)))
    rules=["top500_absfc","top500_p","top500_abst","fc_top500_absfc","p_top500_absfc","pfc_top500_absfc","top250_each","p_top250_each","pfc_top250_each"]
    target=dict(H_only=470,R_only=379,M_only=405,H_R_only=16,R_M_only=16,H_M_only=13,all3=1)
    trials=[]
    for rule in rules:
        selected={s:select(tables[s],rule) for s in CFG}
        mapped={s:(selected[s] if s=="human" else set(orth.loc[orth[s].isin(selected[s]),"human"])) for s in CFG}
        counts=venn(mapped['human'],mapped['rat'],mapped['mouse'])
        trials.append(dict(rule=rule,mapping="frozen one-to-one core; human set not restricted to ortholog universe",selected_counts={s:len(selected[s]) for s in CFG},mapped_counts={s:len(mapped[s]) for s in CFG},ITPR3_selected={s:CFG[s]['gene'] in selected[s] for s in CFG},regions=counts,L1_source_region_distance=sum(abs(counts[k]-target[k]) for k in target),all_three_human_ids=sorted(set.intersection(*mapped.values()))))
    result=dict(scope="Fresh deposit-byte verification and bounded deterministic rerun, not a new DE model, not exact S7A reproduction, not phenotypic/mechanistic efficacy assessment",runtime=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,requests=requests.__version__),source_manifest_sha256=sha(OUT/'source_manifest.json'),source_commit="cdbda819817b7f1005c7b0dcfa2489468347b10e",script_sha256=sha(Path(__file__)),ortholog_map=dict(path=str(MAP.relative_to(ROOT)),sha256=sha(MAP),n=len(orth),itpr3_triple=orth_itpr3.to_dict(orient="records")),diagnostics=diagnostics,itpr3=itpr3,pairwise=pairwise,S7A_source_regions=target,S7A_frozen_map_trials=trials,not_rerun=["ranked GSEA","10,000-draw overlap null","full LOO grid","source limma fit","historical all-homolog live mapping"])
    save_json(OUT / "verification.json",result)
    print(json.dumps(dict(itpr3=itpr3,diagnostics=diagnostics,pairwise=pairwise),indent=2,allow_nan=False))

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--offline",action="store_true")
    parser.add_argument("--output-dir",default="results")
    args=parser.parse_args()
    OFFLINE=args.offline
    OUT=ROOT/args.output_dir
    main()
