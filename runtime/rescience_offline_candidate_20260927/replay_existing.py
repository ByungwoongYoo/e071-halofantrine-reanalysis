"""Replay unmodified upstream diagnostic scripts with frozen local GEO bytes.

The all-homolog BioMart query remains a new dated snapshot, NOT a claim of
byte identity to the unfrozen 2026-09-21 mapping or the paper's gProfiler map.
"""
import contextlib
import argparse
import datetime as dt
import hashlib
import io
import json
from pathlib import Path
import runpy
import requests
import sys

ROOT=Path(__file__).resolve().parent
BASE=ROOT/'upstream/02_COORDINATOR/BIOASSET_DEEP_IP_20260919'
DATA=ROOT/'source_data'
OUT=ROOT/'results'
parser=argparse.ArgumentParser()
parser.add_argument('--offline',action='store_true')
parser.add_argument('--output-dir',default='results')
cli_args=parser.parse_args()
OUT=ROOT/cli_args.output_dir
OUT.mkdir(exist_ok=True)
real_get=requests.get
log=[]

def cached_get(url,*args,**kwargs):
    if url.startswith('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE327nnn/'):
        parts=url.rstrip('/').split('/')
        name=parts[-1] if parts[-1]!='suppl' else parts[-2]+'_suppl_index.html'
        p=DATA/name
        assert p.is_file(),str(p)
        r=requests.Response();r.status_code=200;r._content=p.read_bytes();r.url=url;r.encoding='utf-8'
        return r
    assert url in ['https://www.ensembl.org/biomart/martservice','https://useast.ensembl.org/biomart/martservice'],url
    if cli_args.offline:
        p=DATA/'ensembl_all_homolog_query_result_20260927.tsv'
        r=requests.Response();r.status_code=200;r._content=p.read_bytes();r.url=url;r.encoding='utf-8'
        return r
    kwargs['timeout']=(20,90)
    r=real_get(url,*args,**kwargs)
    if r.ok and r.text.count('\n')>1000:
        p=DATA/'ensembl_all_homolog_query_result_20260927.tsv'
        p.write_bytes(r.content)
        (DATA/'ensembl_all_homolog_query_20260927.xml').write_text(kwargs['params']['query'],encoding='utf-8')
        log.append(dict(endpoint=url,retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),sha256=hashlib.sha256(r.content).hexdigest(),bytes=len(r.content),scope='New current all-homolog snapshot, not source gProfiler map or proven historical-byte match'))
    return r

requests.get=cached_get
for name,stem,marker in [('e071_itpr3_source_diagnostic.py','original_itpr3_replay','E071_ITPR3_SOURCE_DIAGNOSTIC_JSON'),('e071_venn_rule_diagnostic.py','original_nine_rule_replay','E071_VENN_RULE_DIAGNOSTIC_JSON')]:
    buf=io.StringIO()
    with contextlib.redirect_stdout(buf):
        runpy.run_path(str(BASE/name),run_name='__main__')
    body=buf.getvalue().split(marker+'_BEGIN')[1].split(marker+'_END')[0].strip()
    obj=json.loads(body)
    (OUT/(stem+'.json')).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(stem,body,flush=True)
if not cli_args.offline:
    (OUT/'all_homolog_source_receipt.json').write_text(json.dumps(log,indent=2,allow_nan=False)+'\n',encoding='utf-8')
