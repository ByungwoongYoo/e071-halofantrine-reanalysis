"""Compare retained historical manuscript quantities without changing them."""
from pathlib import Path
import csv, hashlib, json, math
R=Path(__file__).resolve().parent
def load(p): return json.loads((R/p).read_text(encoding='utf-8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
c=load('runs/core_seed1/result.json'); c2=load('runs/core_seed77/result.json')
g=load('runs/gsea_current/result.json'); l=load('runs/loo_frozen/result.json')
h=load('historical/frozen_replay_verification.json')
actual={
 'core_common_universe':c['one2one_common_universe'],
 'core_overlap_observed_ge2':c['author_overlap']['observed_ge2'],
 'core_overlap_observed_all3':c['author_overlap']['observed_all3'],
 'core_author_sizes':c['author_overlap']['sizes'],
 'gsea_common_gene_universe':g['common_gene_universe'],
 'gsea_kegg_pathways':g['kegg_pathways'],
 'gsea_FDR05_all3_terms':g['FDR05_all3_terms'],
 'gsea_FDR05_ge2_terms':g['FDR05_ge2_terms']}
for model in c['models']:
 for key in ['tested','nominal_p05_fc12','fdr_selected']:actual[f'core_{model}_{key}']=c['models'][model][key]
 actual[f'gsea_{model}_FDR05_count']=g['models'][model]['FDR05_count']
for pair,x in c['all_gene_pair_metrics'].items():actual[f'core_{pair}_spearman']=x['spearman']
for label,model in [('human','CEP290_LCA'),('p23h','P23H'),('rd10','rd10')]:
 actual[f'gsea_oxphos_{label}_nes']=g['highlighted']['OXIDATIVE PHOSPHORYLATION'][0]['NES_'+model]
checks=[]
for name,x in h['checks'].items():
 expected=x['expected']; value=actual[name]
 tol=x.get('tol',0)
 ok=math.isclose(value,expected,rel_tol=tol,abs_tol=tol) if tol else value==expected
 checks.append(dict(name=name,actual=value,historical=expected,tolerance=tol,match=ok,
                    difference=value-expected if isinstance(value,(float,int)) and isinstance(expected,(float,int)) else None))
rounded=[]
def check_round(name,value,old,decimals):
 tol=0.5*10**(-decimals)+1e-14
 rounded.append(dict(name=name,actual=value,historical_display=old,rounding_tolerance=tol,
                     difference=value-old,match=abs(value-old)<=tol))
for model,values in {'CEP290_LCA':[.9785,.9853,.9568,.9717],
                     'P23H':[.7971,.9657,.8007,.9618],
                     'rd10':[.6788,.9118,.8027,.9045]}.items():
 for key,old in zip(['loo_all_spearman_min','loo_all_spearman_median','loo_high_effect_spearman_min','loo_high_effect_spearman_median'],values):
  check_round(model+'_'+key,l['models'][model][key],old,4)
for pair,values in {'CEP290_LCA__P23H':[-.0429,-.2439,-.0436,.1982],
                    'CEP290_LCA__rd10':[.3116,.0844,.3116,.4221],
                    'P23H__rd10':[.0686,-.1392,.0628,.2844]}.items():
 x=l['cross_model_loo'][pair]
 for key,val,old in zip(['full','min','median','max'],[x['full']['spearman'],x['min'],x['median'],x['max']],values):check_round(pair+'_'+key,val,old,4)
for model,old in [('CEP290_LCA',.01587),('P23H',.19841),('rd10',.3)]:
 check_round(model+'_exact_label_p',l['models'][model]['exact_label_permutation']['exact_p_ge'],old,5)
highlight_checks=[]
with (R/'historical/table3_source_pathway_robustness.csv').open(encoding='utf-8-sig',newline='') as f:
 for row in csv.DictReader(f):
  matches=[x for rows in g['highlighted'].values() for x in rows if x['term']==row['pathway']]
  assert len(matches)==1
  for old_field,display in row.items():
   if old_field=='pathway': continue
   prefix,metric=old_field.split('_'); model={'human':'CEP290_LCA','P23H':'P23H','rd10':'rd10'}[prefix]
   value=matches[0][metric+'_'+model]; old=float(display)
   decimals=len(display.partition('.')[2]); tol=0 if old in (0.,1.) and metric=='FDR' else 0.5*10**(-decimals)+1e-14
   highlight_checks.append(dict(pathway=row['pathway'],field=old_field,current=value,historical_display=old,
     difference=value-old,rounding_tolerance=tol,match=abs(value-old)<=tol))
overlap=[]
for key,old,d in [('null_ge2_mean',56.91,2),('null_ge2_sd',6.86,2),('p_ge2',.630,3),('null_all3_mean',.677,3),('p_all3',.495,3)]:
 value=c['author_overlap'][key]
 overlap.append(dict(field=key,current=value,historical_display=old,difference=value-old,rounds_to_historical=round(value,d)==old))
result=dict(status='PASS_RETAINED_DETERMINISTIC_METRICS' if all(x['match'] for x in checks+rounded+highlight_checks) else 'DIFFERENCES_REQUIRE_REVIEW',
 historical_checks=checks,loo_and_permutation_rounded_checks=rounded,
 highlighted_pathway_display_checks=highlight_checks,
 overlap_historical_display_comparison=overlap,
 hashseed_determinism={'pythonhashseed_values':['1','77'],'full_core_json_equal':c==c2,
 'byte_identical':(R/'runs/core_seed1/result.json').read_bytes()==(R/'runs/core_seed77/result.json').read_bytes(),
 'finding':'Original and frozen overlap_stats already sort the universe. Earlier generic set-order explanation is unsupported for these versions; no scientific formula change was made.'},
 library_provenance=load('inputs/LIBRARY_RECEIPT.json'),
 scope={'loo_deletions':sum(len(x['loo']) for x in l['models'].values()),
 'cross_model_state_pairs':sum(x['all_state_pairs'] for x in l['cross_model_loo'].values()),
 'label_assignments':{m:x['exact_label_permutation']['n_permutations'] for m,x in l['models'].items()},
 'gsea_models':len(g['models']),'gsea_tested_pathways_per_model':g['kegg_pathways'],
 'gsea_permutations_per_model':2000,'overlap_draws_per_test':10000},
 limits=['Comparison is against archived exact central checks and displayed V16/V17 values, not an unavailable complete historical per-pathway table.',
 'Same-name Enrichr snapshot is frozen now; exact historical library bytes remain unverified even if outputs agree.',
 'Current pandas 2.3.3 differs from historically recorded 3.0.6; all runtime versions are recorded, no silent equivalence asserted.',
 'No source limma refit, exact S7A/S7B contract recovery, full manuscript review, public licensing or submission is established.'])
(R/'FULL_REPLAY_COMPARISON.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
print(json.dumps({'status':result['status'],'historical_checks':len(checks),'historical_failures':[x for x in checks if not x['match']],
 'rounded_checks':len(rounded),'rounded_failures':[x for x in rounded if not x['match']],
 'highlighted_checks':len(highlight_checks),'highlighted_failures':[x for x in highlight_checks if not x['match']],
 'overlap':overlap,'hashseed_equal':c==c2,'scope':result['scope']},indent=2))
