"""Independently recompute and print the top three markers per Leiden cluster."""
from pathlib import Path
import scanpy as sc

adata = sc.read_h5ad(Path(__file__).parent / "data" / "pbmc3k.h5ad")
sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", n_genes=3)
result = adata.uns["rank_genes_groups"]
for cluster in result["names"].dtype.names:
    print(f"Cluster {cluster}:")
    for i, gene in enumerate(result["names"][cluster], 1):
        score = result["scores"][cluster][i - 1]
        logfc = result["logfoldchanges"][cluster][i - 1]
        pval = result["pvals_adj"][cluster][i - 1]
        print(f"  {i}. {gene} (logFC={logfc:.3f}, score={score:.3f}, adjusted p={pval:.3e})")
