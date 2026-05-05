# SeuratESFS

An R package that wraps [ESFS (Entropy Sorting Feature Selection)](https://github.com/aradley/ESFS) for Seurat users. Run the complete ESFS workflow directly on Seurat objects, with automatic GPU acceleration via CUDA (NVIDIA) or MLX/Metal (Apple Silicon) where available.

> **Reference:** Radley A. *et al.* (2026) *Entropy Sorting Feature Selection for single-cell RNA sequencing.*
> Preprint: <https://www.biorxiv.org/content/10.64898/2026.01.26.701684v1>

---

## What is ESFS?

ESFS identifies marker genes and performs feature selection in single-cell RNA-seq data using information-theoretic metrics — without requiring prior cell-type labels or dimensionality reduction. It implements three algorithms:

| Algorithm | Purpose |
|-----------|---------|
| **ES-GSS** | Rank genes by pairwise entropy-based correlation |
| **ES-CCF** | Find optimal cluster combinations per gene |
| **ES-FMG** | Select N non-redundant marker genes |

All scientific computation is performed by the [ESFS Python package](https://github.com/aradley/ESFS). This R package is a thin Seurat-idiomatic wrapper using [reticulate](https://rstudio.github.io/reticulate/).

---

## Installation

### 1. Install SeuratESFS

```r
remotes::install_github("aradley/SeuratESFS")
```

### 2. Set up the Python environment (once per machine)

```r
library(SeuratESFS)
install_esfs()   # auto-detects MLX on Apple Silicon, CUDA on NVIDIA, CPU otherwise
```

After installation, add to your `~/.Rprofile`:

```r
options(esfs.conda_env = "r-esfs")
```

### 3. Verify

```r
check_esfs()
```

---

## Quick start

> **Important:** load `SeuratESFS` before Seurat to ensure the Python environment is configured before any other package initialises Python.

```r
library(SeuratESFS)
library(Seurat)

# Scale the count matrix
pbmc3k <- CreateESFSScaledMatrix(pbmc3k)

# Calculate pairwise ES matrices (GPU-accelerated where available)
pbmc3k <- CalcESMatrices(pbmc3k, label = "Self")

# Rank genes by network connectivity
pbmc3k <- ESRankGenes(pbmc3k)

# Visualise top-ranked genes in a UMAP
gene_umap <- PlotESFSGeneUMAP(pbmc3k, top_ranked_genes = 3000, clustering = "hdbscan")

# Compute and plot cell UMAPs per gene cluster
cell_umaps <- GetESFSCellUMAPs(pbmc3k,
                                gene_clust_labels = gene_umap$labels,
                                top_ess_genes     = gene_umap$genes)
PlotESFSCellUMAPs(pbmc3k, cell_umaps[[1]], cell_umaps[[2]], cell_label = "celltype")

# Select 100 non-redundant marker genes
markers <- RunESFMG(pbmc3k, N = 100L)
```

See `vignette("esfs")` or [docs/esfs.Rmd](docs/esfs.Rmd) for the full workflow.

---

## GPU acceleration

| Operation | CUDA (NVIDIA) | MLX (Apple Silicon) | CPU |
|-----------|:---:|:---:|:---:|
| ES matrix calculations | ✓ GPU | ✓ GPU | Numba JIT |
| ES-CCF / ES-FMG | CPU | CPU | CPU |
| Cell UMAP (KNN + layout) | ✓ cuML GPU | CPU | CPU |
| KNN smoothing | ✓ cuML GPU | CPU | CPU |

```r
GetESFSBackendInfo()   # check current backend
UseESFSMLX()           # force Apple Silicon Metal
UseESFSGPU()           # force NVIDIA CUDA
UseESFSCPU()           # force CPU
```

---

## Function reference

| R Function | ESFS Python function | Description |
|------------|---------------------|-------------|
| `CreateESFSScaledMatrix()` | `create_scaled_matrix()` | Scale counts to [0,1] |
| `CalcESMatrices()` | `parallel_calc_es_matrices()` | Pairwise ES-GSS computation |
| `ESRankGenes()` | `ES_rank_genes()` | Rank genes by network weight |
| `PlotESFSGeneUMAP()` | `plot_top_ranked_genes_UMAP()` | Gene-space UMAP |
| `GetESFSCellUMAPs()` | `get_gene_cluster_cell_UMAPs()` | Cell UMAPs per gene cluster |
| `PlotESFSCellUMAPs()` | `plot_gene_cluster_cell_UMAPs()` | Plot cell UMAPs |
| `RunESCCF()` | `ES_CCF()` | Combinatorial cluster finder |
| `RunESFMG()` | `ES_FMG()` | Non-redundant marker gene selection |
| `KNNSmoothESFS()` | `knn_smooth_gene_expression()` | KNN expression smoothing |
| `SaveESFSUMAPModel()` | `save_umap_model()` | Save UMAP model to disk |
| `LoadESFSUMAPModel()` | `load_umap_model()` | Load saved UMAP model |
| `install_esfs()` | — | One-time Python environment setup |
| `check_esfs()` | — | Diagnostic report |

---

## License

GPL-3 © Arthur Radley
