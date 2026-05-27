# scMultiome Pipeline

End-to-end single-cell multiome analysis pipeline for 10X Genomics Multiome data (RNA + ATAC).
Built for HPC environments using LSF (default) or SLURM.

---

## Pipeline Overview

The pipeline is organized into three phases, each submitted as a single HPC job via `run/runmultiome`.

| Phase | Command | Steps | Description |
|-------|---------|-------|-------------|
| 1 | `preprocess` | 1–4 | Init → Create objects → Call peaks → QC |
| 2 | `filter_merge_reduce` | 5–8 | Filter → Merge → Normalize/Reduce → Cluster markers |
| 3 | `annotate_recall_peaks` | 9–10 | Annotate cell types → Recall peaks per cell type |
| 4 | `linkpeaks` | 11 | Link peaks to genes in parallel per group → merge results |

**User checkpoints are required between each phase** — see the walkthrough below.

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `scripts/01_init.R` | Copy raw CellRanger output into `data/raw/` |
| 2 | `scripts/02_create.R` | Create Seurat objects, compute QC metrics (TSS enrichment, nucleosome signal), detect doublets via scDblFinder |
| 3 | `scripts/03_callpeaks.R` | Call peaks with MACS2, replace 10X peak set with custom peaks assay |
| 4 | `scripts/04_qc.R` | Generate QC plots (PDF), doublet summary CSV, write `configs/qc_config.yml` |
| 5 | `scripts/05_filter.R` | Filter cells per sample using thresholds from `qc_config.yml` |
| 6 | `scripts/06_merge.R` | Build consensus peak set, SCTransform per sample, merge all objects |
| 7 | `scripts/07_normalize_reduce.R` | JoinLayers, PCA, TF-IDF/LSI, Harmony (optional), WNN UMAP, multi-resolution clustering |
| 8 | `scripts/08_cluster_markers.R` | Differential expression per cluster via Libra, saves marker tables and annotation plots |
| 9 | `scripts/09_annotate.R` | Left-join cell type labels from `configs/annotations.csv`, remove Discard cells |
| 10 | `scripts/10_recall_peaks.R` | Split by sample × cell type, recall MACS2 peaks per group, build consensus, re-merge |
| 11 | `scripts/11_linkpeaks.R` | Split object by grouping column, submit one LSF job per group to run `LinkPeaks()` in parallel, then submit a merge job to reassemble results onto the full object |

---

## Requirements

- R 4.4.1 (loaded via `module load R/4.4.1`)
- MACS2 or MACS3 (set path in `run/runmultiome`)
- R packages managed via [renv](https://rstudio.github.io/renv/) — no manual installs needed after `init`

---

## Setup (run once per clone)

**1. Clone the repo**
```bash
git clone https://github.com/Alex-Prystupa-Bioinformatics-Projects/scmultiome-pipeline.git
cd scmultiome-pipeline
```

**2. Edit `run/runmultiome` — set your environment variables**
```bash
export project_prefix=myproject          # output file prefix
export my_macs_path=/path/to/macs2       # path to MACS2/3 executable
```

**3. Initialize project directories and restore R environment**
```bash
bash run/runmultiome init
```
Sets up directories and restores the R environment via renv.
Must be run on the **login node** (requires internet access). Takes a few minutes.

**4. Fill in `configs/samplesheet.csv`** — one row per sample:
```
SampleID,path
sample1,/path/to/cellranger/output/sample1
sample2,/path/to/cellranger/output/sample2
```
Additional metadata columns (e.g. `Condition`, `Timepoint`) can be added and will be carried into the Seurat object metadata.

**5. Fill in `configs/pipeline_config.yml`** — set species and genome:
```yaml
# Mouse
species: mouse
genome: mm10
ensdb: EnsDb.Mmusculus.v79
blacklist: blacklist_mm10

# Human (swap to these)
# species: human
# genome: hg38
# ensdb: EnsDb.Hsapiens.v86
# blacklist: blacklist_hg38_unified

# Batch correction — set to a metadata column name, or "none" to skip
harmony_vars: none
```

---

## Running the Pipeline

### Phase 1 — Preprocess (Steps 1–4)

```bash
bash run/runmultiome preprocess
```

Runs: init → create objects → call MACS2 peaks → QC plots

**After Phase 1 — USER ACTION REQUIRED:**
1. Review `output/plots/{prefix}-qc-plots.pdf` — check per-sample QC distributions
2. Review `output/tables/{prefix}-doublet-summary.csv`
3. Edit `configs/qc_config.yml` — set per-sample thresholds (generated automatically, defaults are conservative):
   ```yaml
   sample1:
     min_nCount_RNA: 500
     max_nCount_RNA: 25000
     min_nCount_peaks: 500
     max_nCount_peaks: 100000
     min_TSS: 2
     min_nucleosome: 0
     max_nucleosome: 2
     remove_doublets: false
   ```
4. Optionally set `harmony_vars` in `configs/pipeline_config.yml` if batch correction is needed

---

### Phase 2 — Filter, Merge, Normalize, Reduce (Steps 5–8)

```bash
bash run/runmultiome filter_merge_reduce
```

Runs: filter cells → consensus peaks + merge → normalize/reduce (PCA, LSI, Harmony, WNN UMAP, clustering) → cluster markers

**After Phase 2 — USER ACTION REQUIRED:**
1. Review `output/reports/{prefix}-clustering-report.pptx` — UMAP plots at all resolutions
2. Review `output/markers/{prefix}-08-markers/` — marker gene tables and dot plots per resolution
3. Fill in `configs/annotations.csv` — assign cell type labels to cluster IDs:
   - First column name must match your chosen resolution (e.g. `wsnn_res.0.4`)
   - Add annotation columns with any names you want (e.g. `Annotation_Broad`, `Annotation_Narrow`)
   - Label any clusters to exclude as `Discard`

   Example:
   ```
   wsnn_res.0.4,Annotation_Broad
   0,Fibroblasts
   1,Macrophages
   2,Discard
   3,T_cells
   ```

---

### Phase 3 — Annotate + Recall Peaks (Steps 9–10)

```bash
bash run/runmultiome annotate_recall_peaks
```

Runs: annotate cell types (step 9) → cell-type-aware peak recall (step 10)

Step 10 splits the object by sample × cell type, re-runs MACS2 per group, and builds a new consensus peak set. This enriches the peaks assay with cell-type-specific peaks that would be missed by bulk per-sample calling.

**Output:** `output/RDS-files/{prefix}-10-recall-peaks-obj.RDS` — ready for downstream analyses (link peaks to genes, SCENIC+, etc.)

---

### Phase 4 — Link Peaks to Genes (Step 11)

```bash
bash run/runmultiome linkpeaks
```

Splits the step-10 object by `linkpeaks_grouping_col`, submits one LSF job per group to run Signac `LinkPeaks()` in parallel, then automatically submits a merge job (with LSF dependency) that reassembles all per-group links back onto the full object.

**Before running — configure in `configs/pipeline_config.yml`:**
```yaml
# Split strategy — one LSF job per unique value in this column
linkpeaks_grouping_col: Annotation_Broad   # or: none, orig.ident, Sample_Celltype

# Gene set for LinkPeaks
linkpeaks_genes: all                       # or: variable_features (faster, good for testing)
```

**Job structure:**
- Orchestrator job (lightweight): splits object, submits N child jobs + 1 merge job
- Child jobs (16 cores, 64GB, 36hr each): one per group, saves `output/RDS/linkpeaks/{group}-linkpeaks.RDS`
- Merge job (8 cores, 64GB, 6hr): fires automatically via `#BSUB -w "done(linkpeaks_*)"` — extracts links from all groups, combines into a single GRanges, and assigns back onto the full object

**Monitor:**
```bash
bjobs -J 'linkpeaks_*'
```

**Output:** `output/RDS-files/{prefix}-11-linkpeaks-obj.RDS`

---

## RDS Checkpoints

Each step saves a checkpoint to `output/RDS-files/`. Naming convention:

```
{prefix}-{step_number}-{step_name}-obj[-list].RDS
```

| File | Content |
|------|---------|
| `{prefix}-02-create-obj-list.RDS` | Per-sample Seurat objects, raw |
| `{prefix}-03-callpeaks-obj-list.RDS` | Per-sample objects with MACS2 peaks assay |
| `{prefix}-05-filter-obj-list.RDS` | Per-sample objects after QC filtering |
| `{prefix}-06-merge-obj.RDS` | Single merged object, consensus peaks, SCTransform |
| `{prefix}-07-normalize-reduce-obj.RDS` | Merged object with WNN UMAP + clustering |
| `{prefix}-08-markers-obj.RDS` | Same as step 7 with cluster marker metadata |
| `{prefix}-09-annotate-obj.RDS` | Annotated object, Discard cells removed |
| `{prefix}-10-recall-peaks-obj.RDS` | Annotated object with cell-type-boosted peaks assay |
| `{prefix}-11-linkpeaks-obj.RDS` | Full object with peak-gene links from all groups combined |

---

## Project Structure

```
scmultiome-pipeline/
├── configs/
│   ├── pipeline_config.yml    # species/genome + harmony settings
│   ├── samplesheet.csv        # sample IDs and CellRanger paths
│   ├── qc_config.yml          # QC thresholds (generated after step 4)
│   └── annotations.csv        # cluster → cell type labels (filled in by user after step 8)
├── scripts/
│   ├── 01_init.R              # step scripts
│   ├── ...
│   ├── 10_recall_peaks.R
│   ├── 11_linkpeaks.R         # linkpeaks orchestrator — splits object, submits child + merge jobs
│   ├── functions.R            # core function library
│   ├── genome_utils.R         # species-agnostic genome loader
│   ├── utils.R                # QC yaml writer + filter helpers
│   └── helpers/
│       ├── run_linkpeaks_worker.R  # per-group LSF worker — runs LinkPeaks on one group
│       ├── merge_linkpeaks.R       # merge job — combines per-group links onto full object
│       └── job_merge.sh            # static LSF job script for merge (depends on linkpeaks_*)
├── routes/
│   ├── 01_preprocess.sh           # chains steps 1–4
│   ├── 02_filter_merge_reduce.sh  # chains steps 5–8
│   ├── 03_annotate_recall_peaks.sh  # chains steps 9–10
│   └── 04_linkpeaks.sh            # step 11 orchestration
├── run/
│   └── runmultiome            # HPC job submission entry point
├── output/
│   ├── RDS-files/             # Seurat object checkpoints
│   ├── RDS/
│   │   └── linkpeaks/         # per-group linkpeaks objects (step 11 intermediates)
│   ├── plots/                 # QC and UMAP plots
│   ├── tables/                # QC and doublet summary tables
│   ├── markers/               # cluster marker outputs (step 8)
│   └── reports/               # PPTX clustering report (step 8)
├── renv.lock                  # locked R package versions
└── renv/                      # renv activation files (library is gitignored)
```

---

## Notes

- Seurat objects are named `seu_obj`, `seu_obj_list`, etc. throughout
- Species-agnostic — supports mouse (mm10) and human (hg38) via `configs/pipeline_config.yml`
- HPC scheduler is auto-detected (LSF or SLURM)
- All plots saved as PDF via `pdf()` / `dev.off()`
- Steps 5–6 in Phase 2 are skipped automatically if the merged RDS already exists (saves time on re-runs of step 7)
