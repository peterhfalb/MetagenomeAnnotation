# MetagenomeAnnotation

Functional and taxonomic annotation pipeline for soil metagenome assemblies, focused on eukaryotes (fungi). Follows the NMDC metagenome annotation workflow adapted for SLURM on UMN MSI agate.

Picks up from the output of `../MetagenomeAssembly/`.

---

## Study context

~30–40 soil samples across a pine → oak → oak-maple mesophication gradient (~10 samples per forest type). Annotation targets ecological function: carbon acquisition (CAZymes), N/C cycling pathways (KEGG), pathogenicity (PHI-base), and taxonomic composition (DIAMOND vs NCBI NR).

---

## Prerequisites

**MSI modules used (already available on agate):**

| Module | Used in |
|---|---|
| `metaeuk/6-a5d39d9-gcc-8.2.0-ji6jath` | Step 2 |
| `diamond/2.0.15-gcc-8.2.0-gkldzx7` | Steps 5, 6, database setup |
| `samtools/1.21` | Step 1 |

**Conda environments (created once via `setup_conda_envs.sh`):**

| Environment | Tools |
|---|---|
| `metaG_tiara` | Tiara (PyTorch-based contig classifier — isolated to avoid solver conflicts) |
| `metaG_annotation` | KOfamScan, dbCAN3, featureCounts (subread), MMseqs2, HMMER |

**BLAST+ module** — needed only for the one-time DIAMOND NR database build in step 0. Find the correct module name with `module spider blast` and update `BLAST_MODULE=` in `00_setup_databases.sh` before running.

**R packages** (for step 8):

```r
install.packages(c("data.table", "tidyr", "stringr", "optparse"))
install.packages("taxonomizr")
```

---

## Pipeline steps

### Step 0 — One-time setup (run once, not per project)

```bash
# 1. Create conda environments (login node, interactive, ~10 min)
bash setup_conda_envs.sh

# 2. Download and index all databases (~24–48 hr dominated by DIAMOND NR build)
#    Before submitting: update BLAST_MODULE= in 00_setup_databases.sh
#    PHI-base must be downloaded manually — see Section 5 in 00_setup_databases.sh
bash submit_setup_databases.sh
```

Databases installed to `/projects/standard/kennedyp/shared/databases/metaG_annotation/`:

| Database | Size | Purpose |
|---|---|---|
| DIAMOND NR (built from MSI's BLAST copy) | ~80 GB | Protein-level taxonomy |
| NCBI taxdump + prot.accession2taxid | ~6 GB | Taxonomy tree for LCA |
| KOfam HMM profiles + ko_list | ~5 GB | KEGG KO assignment |
| dbCAN3 HMM + DIAMOND databases | ~1 GB | CAZyme annotation |
| PHI-base DIAMOND database | ~10 MB | Pathogenicity genes |
| Fungi RefSeq proteomes (MMseqs2 DB) | ~5 GB | MetaEuk gene prediction target |

---

### Step 1 — Classify contigs by domain (Tiara)

```bash
bash submit_classify_contigs.sh \
  --assembly-dir  /path/to/MetagenomeAssembly/output \
  --annotation-dir /path/to/annotation/output
```

Each contig is labeled `eukaryota / prokarya / mitochondria / plastid / unknown` by a deep-learning k-mer classifier. Eukaryotic and prokaryotic contigs are split into separate FASTA files for gene prediction.

- Runtime: ~1–4 hr/sample, 16 CPUs, 32 GB
- Outputs: `tiara/`, `euk_contigs/`, `prok_contigs/`

---

### Step 2 — Eukaryotic gene prediction (MetaEuk)

```bash
# Capture the Tiara job ID from step 1 output, then:
bash submit_metaeuk.sh \
  --assembly-dir  /path/to/assembly \
  --annotation-dir /path/to/annotation \
  --after <TIARA_JOB_ID>
```

MetaEuk predicts gene models on eukaryotic contigs via homology to the fungi RefSeq reference database. Handles introns and eukaryotic splice signals that Prodigal cannot.

- Runtime: ~4–12 hr/sample, 32 CPUs, 128 GB
- Outputs per sample: `metaeuk/<SAMPLE>/<SAMPLE>.fas` (proteins), `.gff` (coordinates), `.codon.fas` (CDS nucleotides)

> **Note:** Prokaryotic gene prediction (Prodigal on `prok_contigs/`) is not included in the current pipeline but the classified prokaryotic FASTA files are retained for future use.

---

### Steps 3–7 — Functional annotation and quantification (run in parallel)

All five steps depend only on MetaEuk completing (step 2). Submit them together after step 2 finishes, using `--after <METAEUK_JOB_ID>`.

```bash
METAEUK_JOB_ID=<job id from step 2>
ANN=/path/to/annotation

bash submit_kofam.sh          --annotation-dir $ANN --after $METAEUK_JOB_ID
bash submit_dbcan.sh          --annotation-dir $ANN --after $METAEUK_JOB_ID
bash submit_phibase.sh        --annotation-dir $ANN --after $METAEUK_JOB_ID
bash submit_diamond_nr.sh     --annotation-dir $ANN --after $METAEUK_JOB_ID
bash submit_featurecounts.sh  --assembly-dir /path/to/assembly \
                               --annotation-dir $ANN --after $METAEUK_JOB_ID
```

| Step | Tool | What it produces | Runtime |
|---|---|---|---|
| 3 | KOfamScan | KEGG KO assignments (N/C acquisition pathways) | ~4–8 hr, 32 CPU |
| 4 | dbCAN3 | CAZyme family calls (carbon degradation enzymes) | ~2–4 hr, 16 CPU |
| 5 | DIAMOND vs PHI-base | Pathogenicity gene hits with phenotype labels | ~1 hr, 8 CPU |
| 6 | DIAMOND vs NCBI NR | Per-protein taxonomy (fungi / kingdom / order / genus) | ~8–16 hr, 32 CPU |
| 7 | featureCounts | Raw read counts per predicted gene (uses existing BAM files) | ~15 min, 8 CPU |

**Step 6 note:** DIAMOND vs NR is the slowest step. Concurrency is capped at 3 simultaneous jobs to avoid saturating shared storage I/O on the NR database.

---

### Step 8 — Integration (R)

Once all annotation steps are complete, merge everything into analysis-ready tables:

```bash
Rscript 08_integrate.R \
  --annotation-dir /path/to/annotation \
  --assembly-dir   /path/to/assembly \
  --db-dir         /projects/standard/kennedyp/shared/databases/metaG_annotation
```

Or load interactively by setting variables at the top of the script and sourcing it.

**Outputs in `integrated/`:**

| File | Description |
|---|---|
| `ko_count_matrix.tsv` | KO × sample raw count matrix → DESeq2 / vegan |
| `cazy_count_matrix.tsv` | CAZy family × sample matrix → DESeq2 / vegan |
| `phi_count_matrix.tsv` | PHI-base phenotype × sample matrix |
| `gene_annotations.tsv` | One row per gene: contig coords + KO + CAZy + PHI + LCA taxonomy |
| `gene_counts_raw.tsv` | Gene × sample count matrix (sparse — per-sample assemblies) |
| `summary_stats.tsv` | Per-sample annotation rates (% genes with KO, CAZy, etc.) |
| `integrated_data.RData` | All tables as R objects — `load()` directly into analysis scripts |

**Downstream analysis note:** Pass raw counts directly to DESeq2 (it handles normalization internally). For vegan, normalize with `decostand()` after loading.

---

## Output directory structure

```
annotation_output/
  tiara/            Tiara classification tables + per-sample summary stats
  euk_contigs/      Eukaryotic contig FASTAs (input to MetaEuk)
  prok_contigs/     Prokaryotic contig FASTAs (retained for future use)
  metaeuk/          Gene predictions: .fas (proteins), .gff, .codon.fas
  kofam/            KOfamScan detail-tsv + mapper tables
  dbcan/            dbCAN3 overview.txt per sample
  phibase/          DIAMOND vs PHI-base hit tables
  diamond_nr/       DIAMOND vs NR hit tables (taxonomy, top 10 hits/gene)
  featurecounts/    featureCounts output + .summary QC files
  integrated/       Final merged tables (see step 8 above)
  logs/             SLURM stdout/stderr logs, organized by step
  sample_list.txt   Master sample list (copied from assembly output)
```

---

## Taxonomy note

Taxonomy is assigned at two levels:

1. **Contig level (Tiara, step 1):** domain only — eukaryote vs. prokaryote vs. organelle. Used to route contigs to the correct gene caller.

2. **Protein level (DIAMOND vs NR, step 6):** full lineage via LCA across the top 10 hits within 90% of the best bitscore. Expect reliable resolution to phylum (Basidiomycota/Ascomycota) for most fungal genes, order/family for well-studied groups (Agaricales, Hypocreales), and genus for genes with close reference matches. AM fungi (Glomeromycota) will resolve poorly due to sparse reference genomes.

---

## Rerunning individual samples

All job scripts are idempotent — completed outputs are skipped. To rerun specific samples only:

```bash
bash submit_metaeuk.sh \
  --assembly-dir  /path/to/assembly \
  --annotation-dir /path/to/annotation \
  --samples "SAMPLE_01,SAMPLE_07"
```

To test on 2 samples before a full run:

```bash
bash submit_classify_contigs.sh \
  --assembly-dir  /path/to/assembly \
  --annotation-dir /path/to/annotation \
  --test 2
```
