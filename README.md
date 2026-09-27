# Halofantrine transcriptomic reanalysis

Code accompanying Byungwoong Yoo's reanalysis of the deposited transcriptomic data for human CEP290-LCA retinal organoids, P23H rats and rd10 mice in Kim et al., *Halofantrine protects photoreceptors in multiple models of retinal degeneration* (Research Square v1, [10.21203/rs.3.rs-9511432/v1](https://doi.org/10.21203/rs.3.rs-9511432/v1)).

The analysis tests nine reconstructions of Figure S7A, examines the deposited ITPR3/Itpr3 results, and quantifies cross-model effect concordance, sample-deletion sensitivity, overlap and ranked pathway enrichment. It starts from deposited differential-expression tables and normalized CPM matrices. It does not repeat raw-read processing, refit the source limma model, or reassess photoreceptor protection.

## Installation

The verified runtime is **Windows with Python 3.12.10**. Other systems and Python versions are not supported by this release contract. The driver checks six analysis-library versions before execution; the environment lock records all 22 installed packages from the verified run.

From this repository's root, using Python 3.12.10:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r environment\requirements-lock.txt
.\.venv\Scripts\python.exe -m pip check
```

Use the official Python Package Index. The repository does not bundle an interpreter or package binaries. Installation needs internet access. The `.venv` directory belongs outside `runtime/`, whose file inventory is checked exactly.

## Run

The examples below run from the repository root. Choose an existing parent directory for results and a **new, nonexistent** work-directory name for each run. Do not place the work directory inside `runtime/`.

```powershell
$Python = (Resolve-Path '.\.venv\Scripts\python.exe').Path
$Work = 'D:\E071-runs\run-001'

# Check source integrity and print the plan; no downloads or calculations.
& $Python -B .\runtime\reproduce.py --work-root $Work --python $Python

# Fetch the specified external inputs, then run the full analysis.
& $Python -B .\runtime\reproduce.py --work-root $Work --python $Python --execute --fetch-inputs --allow-mutable
```

`--allow-mutable` permits requests to sources that may change. It never waives the recorded byte counts or SHA-256 checks. A changed or unavailable remote file stops the run; the driver does not substitute a newer input. A successful dry-run does not test installed library versions or input availability.

For a previously acquired input tree, replace `--fetch-inputs --allow-mutable` with `--input-cache 'D:\E071\input-cache'`. That directory must contain the paths specified in `runtime/input_manifest.json`; it is not merely a folder of arbitrarily named downloads. Cache mode makes no download requests. Keep third-party data in a private work directory and follow their source terms.

A failed execution can leave a partial work directory. Keep its logs for diagnosis and use a different new directory when retrying. The scientific stages block outbound Python network operations after acquisition; this is not an operating-system sandbox.

## Inputs and provenance

The input contract specifies 15 assets at 16 target paths. Four Ensembl-derived assets are bundled at five paths. The other inputs are six compressed GEO tables, three GEO supplemental-directory HTML snapshots and two pathway-library representations. The downloader retrieves the raw Enrichr library and derives the normalized representation locally.

Source accessions are [GSE327431](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327431), [GSE327432](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327432) and [GSE327434](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327434), under [GSE327435](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327435). Exact URLs, hashes, sizes and mapping recipes are in the input contract.

The core and symbol-bearing one-to-one mapping snapshots were preserved on 21 September 2026; the all-homologue mapping was recovered on 27 September. The preserved KEGG_2021_Human raw library was downloaded on 27 September. Its identity with the unpreserved historical download is not established. Raw or normalized KEGG GMT files, generated pathway memberships and downloaded GEO differential-expression and CPM files are **not distributed here**. See `THIRD_PARTY_NOTICES.md` for attribution and rights boundaries.

## Verification and outputs

On 27 September 2026, the selected-source package completed all seven stages in a new isolated environment using exact cached inputs: bounded reconstruction; two core runs with Python hash seeds 1 and 77; ranked GSEA; sample deletion and label assignments; result comparison; and five figures. Six bounded outputs, four full-analysis JSON results and five SVGs matched the retained baseline byte-for-byte.

In a separate test on the same date, all eleven external input assets were recovered without the old cache: ten official requests plus one locally normalized GMT. Every size and hash matched. This checks current acquisition, not future URL stability. It is separate from the complete cached-input execution test.

The scientific scope is 25 sample deletions, 257 cross-model state pairs, 252/126/20 treatment-label assignments, 277 pathways per model with 2,000 GSEA permutations, and 10,000 draws for each overlap null. The comparison checks 26 core, 27 deletion/permutation, 18 pathway and five overlap-rounding quantities. These automated comparisons test numerical agreement with the retained reference outputs.

A completed run writes `PORTABLE_EXECUTION_RECEIPT.json` at the work-directory root. Results and comparison records are under `rescience_offline_candidate_20260927/full_replay/`; bounded reconstructions are under `rescience_offline_candidate_20260927/runs/bounded_current/`; the five SVGs are under `figures/`. Stage logs remain in the work directory. The pipeline does not visually inspect its own figures.

Failure to recover the source authors' exact S7A/S7B procedure is distinct from reproducing this specified reanalysis. No treatment efficacy or absence of efficacy follows from a software PASS.

## Author and licence

Byungwoong Yoo, Independent Researcher. [ORCID 0009-0002-1797-3100](https://orcid.org/0009-0002-1797-3100).

The project's analysis and packaging code is released under the MIT licence in `LICENSE`. Third-party datasets, software dependencies, fonts and the source article retain their own terms. The code licence does not grant rights in those materials or license the manuscript. OpenAI Codex assisted with analysis scripts, reproducibility checks, and manuscript drafting and editing.
