# PBMC Cluster Navigator

A single-cell RNA-seq navigator for exploring PBMC clusters with UMAP, gene and QC colouring, marker evidence, and cluster analysis.

## What it does

- Plots all cells on the precomputed UMAP.
- Colours cells by any gene or by `n_genes`, `total_counts`, or `pct_mito`.
- Shows a cluster → top-marker-genes table for all eight Leiden clusters, including scores, log-fold-changes, and adjusted p-values.
- Provides an interactive cluster picker and Run control for exploring the analysis.
- Computes cluster centroid distances and highlights cluster 4's nearest neighbours.

## How to run it

Install the dependencies with pip:

```bash
pip install -r requirements.txt
```

Place `pbmc3k.h5ad` at `data/pbmc3k.h5ad`, which is the path expected by the app. Then run:

```bash
uvicorn app:app --port 8000
```

Open <http://localhost:8000>.

For an optional temporary public HTTPS link, run the app first and then use:

```bash
cloudflared tunnel --url http://localhost:8000
```

A cloudflared tunnel link is temporary and is only live while the machine running the app and tunnel remains running.

## The question

Alex, the PBMC data owner, needed to know whether cluster 4 (163 of 2,700 cells) was a novel population or a known cell type before a Friday PI handoff.

## The answer

Cluster 4 is FCGR3A+ (CD16+) non-classical monocytes — a known, well-documented cell type, not novel. Top markers: FCGR3A, FCER1G, LST1, AIF1, COTL1, IFITM2, IFITM3, FTH1 — a single coherent myeloid signature with strong statistical support (scores 18.7-20.5, p down to 1.6e-89). QC (2.6% mito, healthy range) rules out dying cells; one outlier cell at >20% mito is a single low-quality cell, not a cluster-level signal. Nearest neighboring cluster on the UMAP is cluster 1 (CD14+ classical monocytes) — consistent with known classical/non-classical monocyte subtype biology, not evidence of novelty. Confidence: high.

## Interview

I interviewed Alex, the data owner, before building the app, covering the goal, the available data, done-criteria, and traps such as treating UMAP separation as proof of novelty or overlooking doublets. Alex's quoted QC numbers differed slightly from the app-computed ones, confirming the value of checking the data directly rather than trusting recalled figures.

## AI assistance

Gene-to-cell-type interpretation was done with help from Claude (Anthropic); all data, statistics, and computed numbers come directly from the dataset via this app.

## Live demo

Live demo: <tunnel URL: https://durable-robot-luck-descending.trycloudflare.com, if currently running>
