# MetagenomeAnnotation

Functional and taxonomic annotation pipeline for soil metagenome assemblies, focused on eukaryotes (fungi). Follows the NMDC metagenome annotation workflow adapted for SLURM on UMN MSI.

Picks up from the output of `../MetagenomeAssembly/`.

---

## Study context

40 soil samples from ForestGEO plots across a pine → oak → ash/elm mesophication gradient. Annotation targets ecological function: carbon acquisition (CAZymes, Pfam domains), N/C cycling pathways (KEGG), pathogenicity (PHI-base), and taxonomic composition (MMseqs2 vs. Fungi RefSeq, UniRef90, and Mycocosm/Phytozome). Methodology for the Pfam/Mycocosm-Phytozome/pseudo-genome normalization additions is adapted from *"Potential for functional divergence in ectomycorrhizal fungal communities across a precipitation gradient."*

---

## Prerequisites

**MSI modules used (already available):**

| Module | Used in |
|---|---|
| `metaeuk/6-a5d39d9-gcc-8.2.0-ji6jath` | Step 2 |
| `diamond/2.0.15-gcc-8.2.0-gkldzx7` | Step 5, optional step 6b |
| `samtools/1.21` | Step 1 |

**Conda environments (created once via `setup_conda_envs.sh`):**

| Environment | Tools |
|---|---|
| `metaG_tiara` | Tiara (PyTorch-based contig classifier — isolated to avoid solver conflicts) |
| `metaG_annotation` | KOfamScan, dbCAN3, featureCounts (subread), MMseqs2, HMMER |

**R packages** — install once on MSI before running steps 8–9:

```r
# Step 8 (integration)
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "optparse"))
install.packages("taxonomizr")

# Step 9 (QC report)
install.packages(c("ggplot2", "rmarkdown", "knitr", "kableExtra", "scales"))
```

---

## Pipeline steps

### Step 0 — One-time database setup (run once, not per project)

```bash
# 1. Create conda environments (login node, interactive, ~10 min)
bash setup_conda_envs.sh

# 2. Download and index all databases (~24–48 hr, dominated by MMseqs2 DB build)
bash submit_setup_databases.sh
```

Databases installed to `/projects/standard/kennedyp/shared/databases/metaG_annotation/`:

| Database | Size | Purpose |
|---|---|---|
| Fungi RefSeq proteomes (MMseqs2 DB) | ~5 GB | MetaEuk gene prediction target; MMseqs2 taxonomy |
| KOfam HMM profiles + ko_list | ~5 GB | KEGG KO assignment |
| dbCAN3 HMM + DIAMOND databases | ~1 GB | CAZyme annotation |
| PHI-base DIAMOND database | ~10 MB | Pathogenicity genes |
| NCBI taxdump (names.dmp, nodes.dmp) | ~6 GB | Taxonomy rank expansion via taxonomizr |
| Pfam-A HMM profiles + Pfam-A.hmm.dat | ~100 MB compressed | Pfam domain annotation via hmmscan |
| Mycocosm/Phytozome custom MMseqs2 taxonomy DB *(manual download — see below)* | depends on genome panel | Finer fungal subphylum taxonomy + plant-sequence flagging |
| DIAMOND NR (optional; built from MSI BLAST copy) | ~80 GB | Protein-level taxonomy vs. all NCBI NR |

> **Mycocosm/Phytozome requires a manual download** — JGI genomes need a signed Data Use Agreement and cannot be fetched automatically. Register at the [JGI Genome Portal](https://genome.jgi.doe.gov/portal/), download a panel of fungal (Mycocosm) and plant (Phytozome) protein FASTAs, concatenate them, and supply a genome → NCBI taxid mapping TSV. Full instructions are in `00_setup_databases.sh` Section 8. This database is optional/additive — the pipeline runs fully without it (step 6 UniRef90 taxonomy is unaffected); its absence only means step 6c is skipped.

---

### Steps 1–7 — Submit the pipeline

Edit the two paths at the top of `run_annotation_pipeline.sh`, then run it once from the login node:

```bash
# Edit ASM_DIR and ANN_DIR at the top of the file, then:
bash run_annotation_pipeline.sh

# Test on 2 samples first (recommended before a full run):
bash run_annotation_pipeline.sh --test 2

# Skip database build if already done from a previous project:
bash run_annotation_pipeline.sh --skip-step0
```

The individual `submit_*.sh` scripts still work if you need to rerun a single step.

---

### Step 1 — Classify contigs by domain (Tiara)

Each contig is labeled `eukaryota / prokarya / mitochondria / plastid / unknown` by a deep-learning k-mer classifier. Eukaryotic contigs are separated for gene prediction.

- Runtime: ~1–4 hr/sample, 16 CPUs, 32 GB
- Outputs: `tiara/`, `euk_contigs/`, `prok_contigs/`

---

### Step 2 — Eukaryotic gene prediction (MetaEuk)

MetaEuk v6 predicts gene models on eukaryotic contigs via homology to the Fungi RefSeq reference database. Handles introns and eukaryotic splice signals.

- Runtime: ~4–12 hr/sample, 32 CPUs, 128 GB
- Outputs per sample: `metaeuk/<SAMPLE>/<SAMPLE>.fas` (proteins), `.gff` (coordinates), `.codon.fas` (CDS nucleotides)

> **Gene ID format:** MetaEuk gene IDs are normalized to 4 pipe-delimited fields throughout the pipeline: `targetID|contig|strand|lowerBound`. GFF3 IDs are already in this format; FASTA header IDs (9 fields) are truncated to 4 during integration.

> **Prokaryotic gene prediction** (Prodigal on `prok_contigs/`) is not included in the current pipeline but the classified prokaryotic FASTA files are retained for future use.

---

### Steps 3–7 — Functional annotation and quantification (run in parallel after step 2)

All five steps take the MetaEuk protein FASTA (`metaeuk/<SAMPLE>/<SAMPLE>.fas`) as input and run independently of each other.

---

#### Step 3 — KEGG KO assignment (KOfamScan) — 32 CPU, 32 GB, ~4–8 hr

KOfamScan searches predicted proteins against the **KEGG KOfam database** — ~26,000 HMM profiles, one per KEGG Orthology (KO) number. Each profile is built from curated sequences for that KO and comes with a KO-specific score threshold in `ko_list` that balances precision and recall for that functional group.

**How assignment works:**
1. Each protein is searched against all KOfam profiles using HMMER (`hmmsearch` internally).
2. A hit is **significant** (column 1 = `*`) when its HMMER score exceeds that KO's threshold. Thresholds vary by KO — conserved enzyme families have stricter thresholds than rare or divergent ones.
3. A protein can receive multiple significant KO assignments if it encodes a bifunctional enzyme or a domain shared across KOs.

**Output files per sample:**
- `kofam/<SAMPLE>_kofam.tsv` — full detail-tsv: every gene–KO pair (significant and sub-threshold), with HMMER score and e-value. Useful for re-filtering at different thresholds in R.
- `kofam/<SAMPLE>_kofam_mapper.tsv` — two-column table (`gene_id → KO`) of significant hits only. This is the file read by `08_integrate.R`.

---

#### Step 3b — Pfam domain annotation (HMMER) — 32 CPU, 32 GB, ~4–12 hr

`hmmscan` searches predicted proteins against the **Pfam-A HMM profile database** (~20,000+ domain families), using each family's own curated **gathering threshold** (`--cut_ga`) rather than a single fixed E-value or bitscore cutoff across all families. This is the standard, reproducible approach for Pfam domain annotation. Runs additionally alongside KOfamScan (step 3) — both depend only on MetaEuk.

**How assignment works:**
1. Each protein is searched against the full Pfam-A profile library with `hmmscan --cut_ga --domtblout`.
2. Domain-level output (`--domtblout`, not `--tblout`) is used because a single protein can carry multiple distinct Pfam domains at different coordinate ranges (e.g. a fusion protein).
3. Every reported hit has already passed that family's gathering threshold — unlike KOfamScan's detail-tsv, there is no sub-threshold row to filter out downstream.
4. A gene with multiple distinct domains gets multiple rows in the mapper table, the same multi-hit convention used for KO assignments.

**Output files per sample:**
- `pfam/<SAMPLE>_pfam_domtblout.tsv` — full HMMER3 domtblout (every domain hit, with scores/e-values).
- `pfam/<SAMPLE>_pfam_mapper.tsv` — three-column table (`gene_id, Pfam_name, Pfam_accession`), accession version suffix stripped (e.g. `PF01112.21` → `PF01112`) for stability across Pfam releases. This is the file read by `08_integrate.R`, and the source of the Asparaginase (PF01112) pseudo-genome normalization (see Step 8 and the Integrated output files section below).

---

#### Step 4 — CAZyme annotation (dbCAN3) — 16 CPU, 32 GB, ~2–4 hr

dbCAN3 identifies **carbohydrate-active enzymes (CAZymes)** using two independent evidence streams, run together by `run_dbcan CAZyme_annotation`:

| Tool | Database | How it works |
|---|---|---|
| **HMMER** | dbCAN HMM profiles (one per CAZy family) | Detects conserved domain architecture; good for divergent sequences within a family |
| **DIAMOND** | Characterized CAZy protein sequences | Sequence similarity to biochemically verified CAZymes; good for well-studied families |

**High-confidence call:** a gene is a high-confidence CAZyme when supported by ≥ 2 tools. This threshold is applied in step 8 (`--min-tools`, default 2), not here — all calls are kept in the raw output so the threshold can be adjusted without re-running.

**CAZy class** is extracted from the HMMER family call by stripping the numeric suffix: `GH` (glycoside hydrolases — cellulose, hemicellulose, starch degradation), `AA` (auxiliary activities — oxidative lignin/cellulose degradation), `CE` (carbohydrate esterases — deacetylation of plant polymers), `PL` (polysaccharide lyases — pectin degradation), `CBM` (carbohydrate-binding modules — substrate targeting, not catalytic), `GT` (glycosyltransferases — biosynthesis).

**Output file per sample:** `dbcan/<SAMPLE>/overview.txt` (filename varies by dbCAN version) — one row per gene with the HMMER call, DIAMOND call, and `#ofTools` (number of supporting tools).

---

#### Step 5 — Pathogenicity gene annotation (DIAMOND vs PHI-base) — 8 CPU, 16 GB, ~1 hr

PHI-base (Pathogen-Host Interaction database) is a curated collection of experimentally verified pathogenicity, virulence, and effector genes from fungal and oomycete pathogens. A DIAMOND search finds predicted proteins with homology to these known pathogenicity genes.

**How assignment works:**
- `diamond blastp --sensitive` mode (slower but recovers more divergent homologs), E-value ≤ 1e-5, `--max-target-seqs 1` (best hit per query only). One hit per gene is sufficient here — we want the closest known pathogenicity gene match, not a ranked list.
- The subject title (`stitle`) field encodes the PHI-base annotation as a `#`-delimited string: `UniProtID # PHI:XXXX # gene_name # taxid # pathogen_species # phenotype`. The phenotype (field 6) gives the virulence class: `reduced_virulence`, `loss_of_pathogenicity`, `increased_virulence`, `effector`, or `unaffected`. Parsed in step 8.
- A PHI-base hit does **not** mean the gene is necessarily a pathogenicity factor in your organism — it means the protein is homologous to one. Interpret in the context of the organism's known lifestyle.

**Output file per sample:** `phibase/<SAMPLE>_phibase.tsv` — DIAMOND tabular format (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, stitle).

---

#### Step 6 — Protein-level taxonomy (MMseqs2 easy-taxonomy) — 32 CPU, 128 GB, ~2–6 hr

`mmseqs easy-taxonomy` assigns each predicted protein a taxonomic identity by searching against **UniRef90** and computing a lowest common ancestor (LCA) across the top hits. The full pipeline (database creation → search → LCA → output conversion) runs internally as a single command.

**How assignment works:**
1. **Search:** proteins are searched against UniRef90 at sensitivity `-s 6` (scale 1–7.5; 6 is the recommended metagenomics default). UniRef90 clusters UniProt at 90% identity, providing broad taxonomic coverage with manageable size.
2. **LCA mode 2 (2bLCA):** this is MMseqs2's recommended mode for metagenomics. It combines two signals:
   - the **top-scoring hit** (for specificity — resolves well-matched genes to species/genus)
   - the **LCA across all hits within a bitscore window** (for robustness — prevents a single divergent hit from pulling the assignment up to a wrong taxon)
   
   2bLCA produces assignments that are more resolved than pure LCA while being less susceptible to false specificity than top-hit-only assignment.
3. **Output:** one row per protein. `lca_taxid = 0` = unclassified (no hits above threshold). `--tax-lineage 1` adds the full semicolon-separated lineage path (e.g., `Eukaryota;Fungi;Dikarya;Ascomycota;...`).
4. In step 8, `lca_taxid` values are expanded to standard ranks (superkingdom → species) using **taxonomizr** and the local NCBI taxonomy database.

**Memory note:** the UniRef90 prefilter index requires ~306 GB when fully loaded. `--split-memory-limit 100G` chunks the target database and runs ~3 search passes per sample, making it feasible on 128 GB nodes.

**Output file per sample:** `mmseqs_taxonomy/<SAMPLE>_lca.tsv` — columns: `gene_id`, `lca_taxid`, `lca_rank`, `lca_name`, `lineage`.

---

#### Step 6c — Mycocosm/Phytozome taxonomy, additive (MMseqs2 easy-taxonomy) — 32 CPU, 128 GB, ~8 hr

Runs the same `mmseqs easy-taxonomy` mechanics as step 6 (2bLCA, `--tax-lineage`), but against a **custom MMseqs2 taxonomy database** built from JGI Mycocosm (fungal) and Phytozome (plant) genomes instead of UniRef90. This replicates the taxonomy approach from the source methodology paper for identifying ectomycorrhizal fungi: finer resolution to fungal **subphylum** (Agaricomycotina, Pezizomycotina) and explicit flagging of plant-derived sequences.

**This step is purely additive — no filtering is applied anywhere in the pipeline based on this taxonomy.** It runs alongside step 6 (UniRef90), not in place of it: UniRef90 still provides broad kingdom-level taxonomy/contamination checks; this step adds a second, finer-resolution taxonomic layer. All genes are retained regardless of how they classify here. Optional: if the custom database isn't built (see the manual-download note in the database setup table above), this step is skipped and the rest of the pipeline runs unaffected.

**Output file per sample:** `mycocosm_taxonomy/<SAMPLE>_lca.tsv` — identical 5-column schema to step 6's output (`gene_id`, `lca_taxid`, `lca_rank`, `lca_name`, `lineage`). In `08_integrate.R` and `gene_annotations.tsv`, these columns are prefixed `myco_` to distinguish them from the UniRef90 taxonomy columns.

---

#### Step 6b — Broad protein taxonomy, optional (DIAMOND vs NCBI NR) — 32 CPU, ~8–16 hr

An optional step that searches predicted proteins against all of NCBI NR via DIAMOND. Provides broader taxonomic coverage than the UniRef90 search, at the cost of much higher runtime and storage I/O. Useful for characterizing bacterial/archaeal contamination or non-fungal eukaryotes. Not included in the default pipeline; output in `diamond_nr/`.

---

#### Step 7 — Read quantification (featureCounts) — 8 CPU, 16 GB, ~15 min

featureCounts (from the subread package) intersects each sample's sorted BAM file with MetaEuk's GFF to count how many reads overlap each predicted gene. This produces the raw read counts used in differential abundance analysis. No re-mapping is needed — the BAM files from the assembly step already map reads to all contigs.

**How counting works:**

- **GFF auto-detection:** MetaEuk v6 produces `exon` features with a `Parent=` attribute; older versions produce `CDS` features with an `ID=` attribute. The script detects which is present and sets the feature type (`-t`) and grouping attribute (`-g`) accordingly.
- **Overlap handling (`-O --fraction`):** if a read spans multiple adjacent exon features from the same gene model (common with MetaEuk's tiled exon structure), it is distributed fractionally rather than double-counted.
- **Strandedness (`-s 0`):** unstranded — metagenome DNA reads carry no strand orientation information.
- **Paired-end (`-p --countReadPairs`):** counts read pairs, not individual reads, which better reflects library complexity.

**Why most reads are `Unassigned_NoFeatures`:** the BAM contains reads mapped to all assembled contigs, but MetaEuk only predicted genes on Tiara-classified eukaryotic contigs (≥ 500 bp, typically 5–30% of total contigs). Reads mapping to prokaryotic contigs, short contigs, or unannotated intergenic regions are all `Unassigned_NoFeatures`. A low assigned fraction is expected and is not a problem.

**Output files per sample:**
- `featurecounts/<SAMPLE>_counts.txt` — count table: one row per gene model, columns = Geneid, Chr, Start, End, Strand, Length, counts.
- `featurecounts/<SAMPLE>_counts.txt.summary` — assignment status breakdown (Assigned, Unassigned_NoFeatures, Unassigned_Unmapped, Unassigned_MultiMapping, Unassigned_Ambiguity, etc.).

---

### Step 8 — Integration (R)

Once all annotation steps are complete, merge everything into analysis-ready tables:

```bash
Rscript 08_integrate.R \
  --annotation-dir /projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation \
  --assembly-dir   /path/to/MetagenomeAssembly/output \
  --db-dir         /projects/standard/kennedyp/shared/databases/metaG_annotation
```

Optional flags: `--min-tools 2` (dbCAN confidence threshold), `--top-hits 10`, `--bitscore-frac 0.90`, `--pseudocount 1` (pseudo-genome ALR normalization — see below).

---

### Step 9 — QC report (R Markdown)

Renders a parameterized HTML QC report with 8+ figures covering every pipeline step:

```bash
module load pandoc   # or: conda activate metaG_annotation (if pandoc installed there)

Rscript -e "rmarkdown::render(
  '09_qc_report.Rmd',
  params = list(
    ann_dir  = '/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation'
  ),
  output_file = '/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation/qc_report.html'
)"
```

**Figures generated:**

| Figure | Content |
|---|---|
| 1 | Tiara contig classification (% per domain, stacked bar; all contigs including < 500 bp) |
| 1b | Tiara — classified contigs only (≥ 500 bp), same breakdown re-normalized |
| 2 | MetaEuk predicted protein count per sample |
| 3 | Annotation coverage heatmap (% genes annotated by each tool, grouped by stand type) |
| 4 | CAZyme class breakdown per sample (GH, AA, CE, PL, CBM, GT) |
| 4b | KEGG KO functional categories per sample (downloaded from KEGG REST API; cached after first run) |
| 5 | PHI-base phenotype breakdown per sample |
| 6 | MMseqs2 taxonomy — kingdom/phylum composition per sample |
| 7 | featureCounts read assignment (Assigned, NoFeatures, Unmapped, MultiMapping, Ambiguity) |

All figures are saved as PNGs to `{ann_dir}/qc/figures/`. A per-sample summary table with cell color-coding and `qc_summary.csv` are also written to `{ann_dir}/qc/`.

> **Not yet covered in the QC report:** Pfam family breakdown and Mycocosm/Phytozome taxonomy (plant vs. fungi proportions, fungal subphylum breakdown) figures are planned as a fast-follow, once real data exists from steps 3b/6c to visualize and verify against.

**Stand type metadata** is read from `/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation/stakes.csv`. Samples are joined on stake number (extracted from the sample name pattern `GEO-XX-YY_[Cleaned_]SN`) and grouped/colored by `stand_type` (ash or elm / oak / pine / mixed) in all figures.

---

## Integrated output files

All files written to `{ann_dir}/integrated/` by `08_integrate.R`.

---

### `gene_annotations.tsv`

Master annotation table — **one row per predicted gene**, all annotation layers joined. This is the primary file for custom queries (e.g., "all high-confidence GH18 genes in Ascomycota").

| Column | Description |
|---|---|
| `gene_id` | Normalized MetaEuk gene ID: `targetID\|contig\|strand\|lowerBound` |
| `sample` | Sample name |
| `contig` | Contig the gene was predicted on |
| `start`, `end` | Gene coordinates (bp, 1-based) |
| `strand` | `+` or `-` |
| `length` | Gene length in bp |
| `KO` | KEGG KO ID(s), semicolon-separated if multiple (e.g., `K01083;K07024`) |
| `KO_definitions` | KO description(s), semicolon-separated |
| `CAZy_families` | High-confidence CAZy family call(s), semicolon-separated (e.g., `GH18;CBM50`); NA if gene did not meet the `--min-tools` threshold |
| `CAZy_n_tools` | Number of dbCAN tools supporting the call (max across families) |
| `phi_bitscore` | DIAMOND bitscore for best PHI-base hit |
| `phi_accession` | PHI-base accession (e.g., `PHI:2`) |
| `phi_gene` | Gene name in PHI-base (e.g., `BMP1`) |
| `phi_pathogen` | Pathogen species in PHI-base |
| `phi_phenotype` | Virulence phenotype: `reduced_virulence`, `loss_of_pathogenicity`, `increased_virulence`, `effector`, `unaffected`, or `other` |
| `lca_taxid` | NCBI taxon ID of the MMseqs2 LCA assignment (UniRef90, step 6) |
| `lca_rank` | Taxonomic rank of the LCA (e.g., `species`, `genus`, `phylum`) |
| `lca_name` | Taxon name at the LCA rank |
| `lineage` | Full lineage string from MMseqs2 (UniRef90) |
| `tax_superkingdom` … `tax_species` | Standard rank expansion via taxonomizr (7 columns: superkingdom, phylum, class, order, family, genus, species), UniRef90 source |
| `Pfam_accessions` | Pfam domain accession(s), semicolon-separated if multiple (e.g., `PF00187;PF01112`); bare accessions, version suffix stripped |
| `Pfam_names` | Pfam family name(s), semicolon-separated, matching `Pfam_accessions` order |
| `myco_lca_taxid` | NCBI taxon ID of the Mycocosm/Phytozome LCA assignment (step 6c, additive) |
| `myco_lca_rank` | Taxonomic rank of the Mycocosm/Phytozome LCA |
| `myco_lca_name` | Taxon name at the Mycocosm/Phytozome LCA rank |
| `myco_lineage` | Full lineage string from the Mycocosm/Phytozome taxonomy search |
| `myco_tax_superkingdom` … `myco_tax_species` | Standard rank expansion via taxonomizr, Mycocosm/Phytozome source — only resolves correctly if `taxid_mapping.tsv` used valid NCBI taxids (see Step 6c) |
| `myco_is_plant` | `TRUE` if `myco_lineage` contains "Viridiplantae" — informational flag, not a filter |
| `myco_fungal_subphylum` | Fungal subphylum extracted from `myco_lineage` (e.g., `Agaricomycotina`, `Pezizomycotina`) if present; provisional vocabulary, see note below |

Genes missing an annotation layer have `NA` for that layer's columns. All `myco_*` columns are `NA` for every gene if step 6c's database wasn't built (optional/additive layer).

---

### `gene_counts_raw.tsv`

Raw **gene × sample count matrix** from featureCounts. Rows = `gene_id`, columns = sample names, values = read-pair counts.

> Because each sample has its own assembly, every gene belongs to exactly one sample — the matrix is structurally sparse (each row has counts in one column only). Use the family-level matrices below for cross-sample comparisons.

---

### `ko_count_matrix.tsv`

**KO × sample count matrix** — primary input for KEGG-level differential abundance analysis.

- Rows = KEGG KO ID (e.g., `K01083`)
- Columns = sample names
- Values = summed read-pair counts for all genes in that sample assigned to that KO
- Multi-KO genes (rare) contribute their count to each assigned KO

```r
load("integrated_data.RData")
# DESeq2
dds <- DESeqDataSetFromMatrix(ko_matrix[, -1], colData = meta, design = ~stand_type)
# vegan
ko_hell <- decostand(t(ko_matrix[, -1]), method = "hellinger")
```

---

### `cazy_count_matrix.tsv`

**CAZy family × sample count matrix** — primary input for CAZyme differential abundance analysis.

- Rows = CAZy family ID (e.g., `GH18`, `AA9`, `CBM50`)
- Columns = sample names
- Values = summed read-pair counts for high-confidence (≥ `--min-tools`) CAZyme genes
- Genes with multiple family assignments contribute their count to each family

Useful subsets: filter rows by prefix to isolate GH (cellulose/hemicellulose degradation), AA (oxidative/lignin degradation), PL (pectin), CE (deacetylation), or CBM (binding modules).

---

### `phi_count_matrix.tsv`

**PHI-base phenotype × sample count matrix** — for pathogenicity gene analysis.

- Rows = PHI-base phenotype label
- Columns = sample names
- Values = summed read-pair counts for genes with a PHI-base hit in that phenotype category

---

### `pfam_count_matrix.tsv`

**Pfam family × sample count matrix** — built the same way as the KO/CAZy/PHI matrices above, via the same generic `build_family_matrix()` helper.

- Rows = Pfam accession (e.g., `PF01112`)
- Columns = sample names
- Values = summed read-pair counts for genes with that Pfam domain
- Genes with multiple distinct Pfam domains contribute their count to each family

---

### Pseudo-genome (ALR) normalized matrices

`ko_count_matrix_normalized.tsv`, `cazy_count_matrix_normalized.tsv`, `phi_count_matrix_normalized.tsv`, `pfam_count_matrix_normalized.tsv`, and `asparaginase_pseudogenome_counts.tsv`.

This replicates the normalization approach from the source methodology paper: **Asparaginase (Pfam `PF01112`)** is treated as a near-single-copy marker gene, used as a proxy for the genome equivalents represented by each sample's sequencing depth. For each sample, the Asparaginase pseudo-genome count is the summed featureCounts read-pair count across all genes with a PF01112 hit (written out directly in `asparaginase_pseudogenome_counts.tsv` for auditing).

Each of the four family matrices is normalized per sample as:

```
ratio     = (family_count + pseudocount) / (asparaginase_count + pseudocount)
normalized = log(ratio)
```

This is mathematically an **additive log-ratio (ALR) transform** using Asparaginase as the reference component. `pseudocount` (default 1, set via `--pseudocount`) is added to both numerator and denominator to avoid `log(0)`; the "right" value is somewhat arbitrary and depth-sensitive — treat the default as a starting point requiring a sensitivity check (e.g. rerun with `--pseudocount 0.5` or `5` and confirm conclusions are stable), not a fixed constant.

**Why all four matrices, not just Pfam:** the KO, CAZy, PHI, and Pfam matrices are all built from the same underlying featureCounts read-pair counts, just bucketed by different functional-family assignments from different annotation tools — so the numerator and denominator are always on identical units regardless of which tool produced the family label. The Asparaginase count functions as a **sample-level scaling factor** (how much genomic content that sample's sequencing represents), a property of the sample rather than of the annotation system used to derive the numerator — the same logic behind single-copy-marker "genome equivalents" normalization methods (e.g. MicrobeCensus). **This is a methodological extension beyond the cited paper's literal scope**: the paper only had Pfam-derived functional counts to normalize (Pfam was its only functional annotation system), so it never tested this logic against KO/CAZy/PHI-style counts from other annotation tools. Treat the Pfam-normalized matrix as the direct replication and the KO/CAZy/PHI-normalized matrices as an analogous extension.

**Caveats:**
- Samples with a near-zero Asparaginase pseudo-genome count have unreliable normalized values for all four matrices — their ALR values become dominated by the pseudocount rather than reflecting real abundance. `08_integrate.R` prints a warning listing any such samples; check `asparaginase_pseudogenome_counts.tsv` before trusting normalized values for those samples.
- The PF01112 row in `pfam_count_matrix_normalized.tsv` is exactly `log(1) = 0` for every sample, by construction (it's normalized against itself). **This is an expected transform artifact, not a biological signal** — do not interpret it as "Asparaginase abundance is flat across samples."
- The Asparaginase single-copy-marker assumption is inherited from the source paper's methodology and cannot be independently validated by this pipeline.
- The `*_normalized.tsv` files are already log-transformed — do not feed them into DESeq2 or `vegan::decostand()`, which expect raw counts. Use the raw `*_count_matrix.tsv` files for those tools instead.

---

### `summary_stats.tsv`

Per-sample QC metrics table. One row per sample.

| Column | Description |
|---|---|
| `sample` | Sample name |
| `genes_predicted` | Total genes predicted by MetaEuk |
| `total_read_count` | Total read pairs assigned by featureCounts |
| `genes_with_KO` | Genes with at least one KOfamScan hit |
| `pct_KO` | % genes with KO annotation |
| `genes_with_CAZy` | Genes with a high-confidence dbCAN call |
| `pct_CAZy` | % genes with CAZy annotation |
| `genes_with_PHI_hit` | Genes with a PHI-base DIAMOND hit |
| `pct_PHI` | % genes with PHI annotation |
| `genes_with_taxonomy` | Genes with a classified MMseqs2 LCA taxid (UniRef90, taxid ≠ 0) |
| `pct_taxonomy` | % genes with UniRef90 taxonomy assignment |
| `genes_with_Pfam` | Genes with at least one Pfam domain hit |
| `pct_Pfam` | % genes with Pfam annotation |
| `genes_with_myco_taxonomy` | Genes with a classified Mycocosm/Phytozome LCA taxid (additive layer) |
| `pct_myco_taxonomy` | % genes with Mycocosm/Phytozome taxonomy assignment |
| `asparaginase_pseudogenome_count` | Raw PF01112 read-pair count for this sample — the ALR normalization denominator (see above) |

---

### `integrated_data.RData`

All tables above saved as named R objects. Load directly into analysis scripts:

```r
load("/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation/integrated/integrated_data.RData")
# Objects: base_dt, count_matrix, ko_matrix, cazy_matrix, phi_matrix, pfam_matrix,
#          ko_matrix_normalized, cazy_matrix_normalized, phi_matrix_normalized,
#          pfam_matrix_normalized, asparaginase_dt, summary_dt
```

---

## Output directory structure

```
MetaG_Annotation/
  tiara/               Tiara classification tables + per-sample summary stats
  euk_contigs/         Eukaryotic contig FASTAs (input to MetaEuk)
  prok_contigs/        Prokaryotic contig FASTAs (retained for future use)
  metaeuk/             Gene predictions: .fas (proteins), .gff, .codon.fas
  kofam/               KOfamScan mapper TSVs (_kofam_mapper.tsv per sample)
  pfam/                Pfam hmmscan domtblout + mapper TSVs per sample
  dbcan/               dbCAN3 overview file per sample
  phibase/             DIAMOND vs PHI-base hit tables
  mmseqs_taxonomy/     MMseqs2 LCA taxonomy tables vs UniRef90 (_lca.tsv per sample)
  mycocosm_taxonomy/   MMseqs2 LCA taxonomy tables vs Mycocosm/Phytozome (additive, optional)
  diamond_nr/          DIAMOND vs NCBI NR hit tables (optional step 6b)
  featurecounts/       featureCounts count tables + .summary QC files
  integrated/          Final merged tables (see above)
  qc/                  QC report HTML + qc_summary.csv + figures/
  logs/                SLURM stdout/stderr logs, organized by step
  sample_list.txt      Master sample list
  stakes.csv           ForestGEO plot metadata (stake, stand_type, metagenomics)
```

---

## Taxonomy note

Taxonomy is assigned at three levels:

1. **Contig level (Tiara, step 1):** domain only — eukaryote vs. prokaryote vs. organelle. Used to route contigs to the correct gene caller. Tiara requires contigs ≥ 500 bp; shorter contigs are excluded from classification.

2. **Protein level, broad (MMseqs2 vs UniRef90, step 6):** LCA computed natively by MMseqs2 across top hits within 90% of the best bitscore against UniRef90. Expanded to standard ranks (superkingdom → species) via taxonomizr in step 8. Expect reliable resolution to phylum (Ascomycota/Basidiomycota) for most fungal genes, order/family for well-studied groups (Agaricales, Hypocreales), and genus for genes with close reference matches. AM fungi (Glomeromycota) resolve poorly due to sparse reference genomes.

3. **Protein level, fine fungal/plant resolution (MMseqs2 vs Mycocosm/Phytozome, step 6c, additive):** same LCA mechanics as (2), but against a custom database of curated JGI fungal and plant genomes. Purpose-built for two things UniRef90 resolves poorly: fungal **subphylum**-level placement (Agaricomycotina, Pezizomycotina) and explicit flagging of plant-derived sequences (`myco_is_plant`). Runs alongside (2), not in place of it — no filtering is applied based on this layer anywhere in the pipeline; both taxonomy layers are retained side by side in `gene_annotations.tsv` for every gene. Optional: skipped entirely if the custom database isn't built (see Step 6c).

---

## Rerunning individual samples

All job scripts are idempotent — completed outputs are skipped. To rerun specific samples only:

```bash
bash submit_metaeuk.sh \
  --assembly-dir   /path/to/assembly \
  --annotation-dir /path/to/annotation \
  --samples "SAMPLE_01,SAMPLE_07"
```
