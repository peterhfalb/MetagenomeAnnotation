# MetagenomeAnnotation

Functional and taxonomic annotation pipeline for soil metagenome assemblies, focused on eukaryotes (fungi). Follows the NMDC metagenome annotation workflow adapted for SLURM on UMN MSI.

Picks up from the output of `../MetagenomeAssembly/`.

---

## Study context

40 soil samples from ForestGEO plots across a pine → oak → ash/elm mesophication gradient. Annotation targets ecological function: carbon acquisition (CAZymes), N/C cycling pathways (KEGG), pathogenicity (PHI-base), and taxonomic composition (MMseqs2 vs. Fungi RefSeq).

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
| DIAMOND NR (optional; built from MSI BLAST copy) | ~80 GB | Protein-level taxonomy vs. all NCBI NR |

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

Optional flags: `--min-tools 2` (dbCAN confidence threshold), `--top-hits 10`, `--bitscore-frac 0.90`.

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
| `lca_taxid` | NCBI taxon ID of the MMseqs2 LCA assignment |
| `lca_rank` | Taxonomic rank of the LCA (e.g., `species`, `genus`, `phylum`) |
| `lca_name` | Taxon name at the LCA rank |
| `lineage` | Full lineage string from MMseqs2 |
| `tax_superkingdom` … `tax_species` | Standard rank expansion via taxonomizr (7 columns: superkingdom, phylum, class, order, family, genus, species) |

Genes missing an annotation layer have `NA` for that layer's columns.

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
| `genes_with_taxonomy` | Genes with a classified MMseqs2 LCA taxid (taxid ≠ 0) |
| `pct_taxonomy` | % genes with taxonomy assignment |

---

### `integrated_data.RData`

All tables above saved as named R objects. Load directly into analysis scripts:

```r
load("/projects/standard/kennedyp/shared/projects/ForestGEO/MetaG_Annotation/integrated/integrated_data.RData")
# Objects: base_dt, count_matrix, ko_matrix, cazy_matrix, phi_matrix, summary_dt
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
  dbcan/               dbCAN3 overview file per sample
  phibase/             DIAMOND vs PHI-base hit tables
  mmseqs_taxonomy/     MMseqs2 LCA taxonomy tables (_lca.tsv per sample)
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

Taxonomy is assigned at two levels:

1. **Contig level (Tiara, step 1):** domain only — eukaryote vs. prokaryote vs. organelle. Used to route contigs to the correct gene caller. Tiara requires contigs ≥ 500 bp; shorter contigs are excluded from classification.

2. **Protein level (MMseqs2, step 6):** LCA computed natively by MMseqs2 across top hits within 90% of the best bitscore against the Fungi RefSeq proteome database. Expanded to standard ranks (superkingdom → species) via taxonomizr in step 8. Expect reliable resolution to phylum (Ascomycota/Basidiomycota) for most fungal genes, order/family for well-studied groups (Agaricales, Hypocreales), and genus for genes with close reference matches. AM fungi (Glomeromycota) resolve poorly due to sparse reference genomes.

---

## Rerunning individual samples

All job scripts are idempotent — completed outputs are skipped. To rerun specific samples only:

```bash
bash submit_metaeuk.sh \
  --assembly-dir   /path/to/assembly \
  --annotation-dir /path/to/annotation \
  --samples "SAMPLE_01,SAMPLE_07"
```
