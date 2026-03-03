# scMultiome Pipeline

End-to-end single-cell multiome analysis pipeline for 10X Genomics Multiome data (RNA + ATAC).
Built for HPC environments using LSF (default) or SLURM.

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `scripts/01_init.R` | Copy raw CellRanger output into `data/raw/` |
| 2 | `scripts/02_create.R` | Create Seurat objects, compute QC metrics, detect doublets |
| 3 | `scripts/03_callpeaks.R` | Call peaks with MACS2, create peaks assay |
| 4 | `scripts/04_qc.R` | Generate QC plots, doublet summary, write `qc_config.yml` |
| 5 | `scripts/05_filter.R` | Filter cells based on QC thresholds *(in progress)* |

---

## Requirements

- R 4.4.1 (loaded via `module load R/4.4.1`)
- MACS2/3 (loaded via module or conda)
- R packages managed via [renv](https://rstudio.github.io/renv/) — no manual installs needed

---

## Setup (run once per clone)

**1. Clone the repo**
```bash
git clone https://github.com/Alex-Prystupa-Bioinformatics-Projects/scmultiome-pipeline.git
cd scmultiome-pipeline
```

**2. Initialize project directories**
```bash
run/runmultiome init
```

**3. Restore the R environment** *(must be run on the login node — requires internet)*
```bash
module load R/4.4.1 && Rscript -e "renv::install('renv'); renv::restore()"
```
This installs all required R packages into `renv/library/` from the locked versions in `renv.lock`.
You only need to do this once per clone. The library lives on the shared filesystem and is accessible by compute nodes.

**4. Configure your project**

Edit `configs/samplesheet.csv` — one row per sample:
```
SampleID,path
sample1,/path/to/cellranger/output/sample1
sample2,/path/to/cellranger/output/sample2
```

Edit `configs/pipeline_config.yml` — set species and genome:
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
```

---

## Running the Pipeline

**Phase 1 — Preprocessing (Steps 1-4)**
```bash
run/runmultiome preprocess
```
Submits an HPC job. Monitor logs in `scripts/outs/`.

After completion:
1. Review `output/plots/{project_prefix}-qc-plots.pdf`
2. Review `output/tables/{project_prefix}-doublet-summary.csv`
3. Edit `configs/qc_config.yml` with per-sample QC thresholds
4. Run Phase 2 (filter) — *coming soon*

---

## Project Structure

```
scmultiome-pipeline/
├── configs/
│   ├── pipeline_config.yml   # species/genome settings
│   ├── samplesheet.csv       # sample IDs and CellRanger paths
│   └── qc_config.yml         # QC thresholds (generated after step 4)
├── scripts/                  # R scripts for each pipeline step
├── routes/                   # shell scripts that chain steps together
├── run/
│   └── runmultiome           # HPC job submission entry point
├── output/
│   ├── RDS-files/            # Seurat object checkpoints
│   ├── plots/                # QC plots
│   └── tables/               # summary tables
├── renv.lock                 # locked R package versions
└── renv/                     # renv activation files (library is gitignored)
```

---

## Notes

- Seurat objects are named `seu_obj`, `seu_obj_list`, etc. throughout
- RDS files follow the naming convention: `{prefix}-{step_number}-{step_name}-obj-list.RDS`
- The pipeline is species-agnostic — supports mouse (mm10) and human (hg38)
- HPC scheduler is auto-detected (LSF or SLURM)
