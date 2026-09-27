# Third-party inputs and dependencies

The MIT licence in this repository applies to the project-authored analysis and packaging code. It does not replace the terms of the resources below. No third-party software binaries or source-article PDF are distributed in the runtime package.

## Ensembl-derived mapping snapshots

The bundled tables contain project-filtered human, mouse and rat orthology records obtained from Ensembl BioMart. One table also contains external gene symbols. Attribute Ensembl and EMBL-EBI when using these inputs: Yates AD et al., *Ensembl 2026*, [10.1093/nar/gkaf1239](https://doi.org/10.1093/nar/gkaf1239).

The one-to-one ID table contains 15,956 rows, and the symbol-bearing table 15,940. Both were preserved on 21 September 2026. The all-homologue response and query XML were preserved on 27 September. The input contract records exact hashes, source endpoints and filtering rules. The historical live-service download was not hashed at the time; its precise release provenance and byte identity cannot be inferred from equal row counts.

[Ensembl's data disclaimer](https://jun2026.archive.ensembl.org/info/about/legal/disclaimer.html) permits use of project-generated data without restriction, subject to third-party constraints. [EMBL-EBI's terms](https://www.ebi.ac.uk/about/terms-of-use/) address redistribution and attribution. External gene names may originate from other resources; the complete symbol-bearing table is not asserted to be CC0 or covered by this project's MIT licence. Neither organisation is claimed to endorse this reanalysis.

## GEO

The downloaded differential-expression and normalized-expression files belong to GSE327431, GSE327432 and GSE327434. Cite the source study and GEO when using them. They are obtained from official NCBI URLs, not redistributed in this repository. NCBI's [disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) does not grant third-party intellectual-property rights held by data submitters.

## Enrichr and KEGG

The pathway analysis obtains KEGG_2021_Human from the official Enrichr endpoint and normalizes it locally using the recorded recipe. This repository does not distribute raw or normalized GMT files or generated pathway memberships. Access does not imply permission to redistribute or use a resource for every purpose; users should consult the [Enrichr resource](https://maayanlab.cloud/Enrichr/) and [KEGG terms](https://www.kegg.jp/kegg/legal.html).

Citations: Xie Z et al., *Gene set knowledge discovery with Enrichr*, [10.1002/cpz1.90](https://doi.org/10.1002/cpz1.90); Kanehisa M et al., *KEGG for taxonomy-based analysis of pathways and genomes*, [10.1093/nar/gkac963](https://doi.org/10.1093/nar/gkac963).

## Python dependencies

The environment lock lists separately installed dependencies; their licences remain with their respective authors. GSEApy is distributed under BSD-3-Clause, Copyright 2016–2020 Zhuoqing Fang. Its installed distribution carries that licence. GSEApy source or binaries are not bundled here. Cite Fang Z, Liu X and Peltz G, *GSEApy: a comprehensive package for performing gene set enrichment analysis in Python*, [10.1093/bioinformatics/btac757](https://doi.org/10.1093/bioinformatics/btac757).
