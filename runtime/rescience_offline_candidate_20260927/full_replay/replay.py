"""Preserved E071 core/GSEA/LOO replay with explicit frozen input adapters."""
from pathlib import Path
import argparse, ast, contextlib, datetime, hashlib, io, json, os, runpy, shutil, sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PACK = ROOT.parent
BASE = Path('02_COORDINATOR/BIOASSET_DEEP_IP_20260919')
os.environ['MPLCONFIGDIR'] = str(ROOT / 'runtime_cache/matplotlib')
sys.path.insert(0, str(ROOT / 'deps'))
import requests

SCRIPTS = {'core': ('e071_gse_replay_frozen.py', 'E071_RESULT_JSON'),
           'gsea': ('e071_ranked_gsea_frozen.py', 'E071_GSEA_JSON'),
           'loo': ('e071_sample_level_robustness.py', 'E071_SAMPLE_LEVEL_ROBUSTNESS_JSON')}
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p, obj): p.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def prepare():
    original = PACK.parent / 'upstream' / BASE
    mappings = BASE / 'SREP_SUBMISSION_PACKAGE_V1/ortholog_mapping'
    pairs = [(original / ('SREP_REVIEW_BUNDLE_V1/code/' + name), ROOT / 'source' / name)
             for name, _ in SCRIPTS.values()]
    pairs += [(original / 'SREP_REVIEW_BUNDLE_V1/ortholog_mapping' / name,
               ROOT / mappings / name)
              for name in ['ensembl_one2one_core.tsv', 'ensembl_one2one_gsea.tsv']]
    pairs += [(original / 'SREP_REVIEW_BUNDLE_V1/package/frozen_replay_verification.json',
               ROOT / 'historical/frozen_replay_verification.json')]
    pairs += [(p, ROOT/'historical'/p.name) for p in original.glob('1[67]_E071*.md')]
    pairs += [(original/'29_E071_SREP_REFERENCE_REPRO_DRAFT_V29.md', ROOT/'historical/V29.md')]
    record=[]
    for src,dst in pairs:
        assert not dst.exists(), dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src,dst)
        assert sha(src)==sha(dst)
        record.append({'file':dst.relative_to(ROOT).as_posix(),'sha256':sha(dst),'bytes':dst.stat().st_size})
    put(ROOT/'COPIED_INPUTS.json',record)

def acquire():
    dst=ROOT/'inputs'; dst.mkdir(exist_ok=True)
    url='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2021_Human'
    raw=dst/'KEGG_2021_Human.download_20260927.raw.gmt'
    assert not raw.exists()
    response=requests.get(url, timeout=(15,60)); response.raise_for_status()
    text=response.content.decode('utf-8')
    assert '<html' not in text[:500].lower()
    rows=[]; terms=set()
    for line in text.splitlines():
        fields=line.strip().split('\t')
        if len(fields)<2: continue
        assert fields[0] not in terms, fields[0]
        terms.add(fields[0])
        genes=[g.split(',')[0] for g in fields[2:] if g.split(',')[0]]
        rows.append('\t'.join([fields[0],fields[1]]+genes))
    assert len(rows)>250
    raw.write_bytes(response.content)
    normalized=dst/'KEGG_2021_Human.gseapy_normalized.gmt'
    normalized.write_text('\n'.join(rows)+'\n',encoding='utf-8')
    put(dst/'LIBRARY_RECEIPT.json',{'library_name':'KEGG_2021_Human','retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'url':url,'final_url':response.url,'status_code':response.status_code,
        'raw_sha256':sha(raw),'raw_bytes':raw.stat().st_size,'term_count':len(rows),
        'normalized_sha256':sha(normalized),'normalization':'Exact GSEApy 1.3.1 download_libraries parsing: strip lines, split tab, gene text before comma, remove empty genes, preserve order',
        'headers':{k:v for k,v in response.headers.items() if k.lower() in ['date','etag','last-modified','content-type','content-length']},
        'official_endpoint_documentation':'https://github.com/MaayanLab/backr/blob/main/README.md',
        'historical_byte_identity':'NOT_ESTABLISHED: no historical raw GMT or hash was retained; same name and matching outputs do not prove identical bytes',
        'rights_boundary':'Public official download for private reproducibility analysis. No public redistribution or blanket relicensing is performed.'})
    print(json.dumps(json.loads((dst/'LIBRARY_RECEIPT.json').read_text()),indent=2))

def execute(component, run_name):
    import numpy as np, pandas as pd, scipy, socket
    from importlib.metadata import version
    out=ROOT/'runs'/run_name
    assert out.parent.resolve()==(ROOT/'runs').resolve()
    out.mkdir(parents=True,exist_ok=False)
    original_get=requests.get
    def local_get(url,*args,**kwargs):
        assert url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/'), url
        parts=url.rstrip('/').split('/')
        name=parts[-1] if parts[-1]!='suppl' else parts[-2]+'_suppl_index.html'
        path=PACK/'source_data'/name
        assert path.is_file(), path
        r=requests.Response();r.status_code=200;r._content=path.read_bytes();r.url=url;r.encoding='utf-8'
        return r
    requests.get=local_get
    blocked=[]
    def guard(event,args):
        if event=='socket.bind' and args[1] in [('::1',0),('127.0.0.1',0)]: return
        if event in ['socket.connect','socket.getaddrinfo','socket.sendto','socket.bind']:
            blocked.append(event);raise RuntimeError('No analysis network '+event)
    sys.addaudithook(guard)
    if component=='gsea':
        import gseapy as gp
        reader=gp.get_library
        def local_library(name,organism='Human',**kwargs):
            assert name=='KEGG_2021_Human' and organism=='Human'
            return reader(str(ROOT/'inputs/KEGG_2021_Human.gseapy_normalized.gmt'),organism=organism,**kwargs)
        gp.get_library=local_library
    name,marker=SCRIPTS[component]
    path=ROOT/'source'/name
    source=path.read_text(encoding='utf-8')
    tree=ast.parse(source, filename=str(path))
    adapters=['requests.get resolves only hash-preserved local GEO files']
    if component=='loo':
        replacement=ast.parse("def biomart():\n return pd.read_csv(FROZEN_CORE_MAP, sep='\\t', dtype=str).fillna('')\n").body[0]
        found=sum(isinstance(n,ast.FunctionDef) and n.name=='biomart' for n in tree.body)
        assert found==1
        tree.body=[replacement if isinstance(n,ast.FunctionDef) and n.name=='biomart' else n for n in tree.body]
        ast.fix_missing_locations(tree)
        adapters.append('biomart acquisition function reads preserved 15956-row one-to-one map; no arithmetic/selection change')
    if component=='gsea': adapters.append('gp.get_library reads dated official Enrichr snapshot with equivalent GSEApy parsing; unknown historical byte identity')
    capture=io.StringIO(); oldcwd=Path.cwd(); os.chdir(ROOT)
    try:
        ns={'__name__':'__main__','__file__':str(path),'FROZEN_CORE_MAP':str(ROOT/BASE/'SREP_SUBMISSION_PACKAGE_V1/ortholog_mapping/ensembl_one2one_core.tsv')}
        with contextlib.redirect_stdout(capture): exec(compile(tree,str(path),'exec'),ns)
    finally: os.chdir(oldcwd)
    log=capture.getvalue(); (out/'stdout.txt').write_text(log,encoding='utf-8')
    obj=json.loads(log.split(marker+'_BEGIN')[1].split(marker+'_END')[0].strip())
    assert 'biomart_error' not in obj, obj.get('biomart_error')
    assert not blocked, blocked
    put(out/'result.json',obj)
    if component=='gsea':
        for model,df in ns['results'].items(): df.to_csv(out/(model+'.csv'),index=False)
        put(out/'tested_gene_sets.json',ns['lib2'])
    put(out/'RUN_RECEIPT.json',{'component':component,'status':'EXECUTED_OFFLINE','source_sha256':sha(path),'result_sha256':sha(out/'result.json'),
        'runtime':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'gseapy':version('gseapy') if component=='gsea' else None},
        'pythonhashseed':os.environ.get('PYTHONHASHSEED'),'adapters':adapters,'analysis_network_attempts':len(blocked),'calculation_rules_modified':False})
    print(component,sha(out/'result.json'),json.dumps(obj)[:300])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','acquire','core','gsea','loo']);p.add_argument('--run-name',default='run1');a=p.parse_args()
    if a.action=='prepare':prepare()
    elif a.action=='acquire':acquire()
    else:execute(a.action,a.run_name)
