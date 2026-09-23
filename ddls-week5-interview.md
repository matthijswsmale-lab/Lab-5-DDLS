# DDLS 2026 — Week 5 interview transcript

**Lab:** Interview a Single-Cell Data Owner
**Data owner:** Alex Rios — PhD student who thinks they've found a novel immune cell population
**Student:** Matthijs Wessel Smale <smale@kth.se>
**Exported:** 2026-09-23T11:36:45+00:00

> This is the record of your interview with the data owner. Hand it to your analyst agent as the raw source for `AGENTS.md` and `spec.md`.

---

**You**:

hey alex what are you trying to figure out and what happens once you know?

**Alex Rios**:

I’m Alex, and this is my first PBMC single-cell dataset—2,700 cells split into numbered clusters 0–7. **Cluster 4**, about 163 cells, sits apart as its own group and I can’t place it; I’m wondering if it might be a novel population. I need you to help determine whether it’s genuinely new rather than just an unlabelled known type. If it holds up, I’d take it to my PI and plan follow-up validation; if not, I’ll label it properly and move on.

**You**:

what doe the files look like? what is one row of data?

**Alex Rios**:

I don’t have the exact file layout in front of me—just the cluster summary and expression data. Conceptually, one row is one cell: its gene-expression values, plus metadata like cluster assignment, detected genes, and mitochondrial percentage. Do you mean the expression matrix or the cluster-summary table?

**You**:

does cluster four have a cluster name

**Alex Rios**:

No—it’s only labeled **cluster 4** right now. I haven’t assigned it a cell type, which is why I’m hoping the marker pattern can tell us whether it’s **a new population** or something known.

**You**:

Would it satisfy you if I told you cluster 4 is actually a known cell type I can name — or are you only interested if the answer is 'novel'?

**Alex Rios**:

A defensible known identity would satisfy me—novelty isn’t the only acceptable outcome. **Separate doesn’t mean new**, does it? I got ahead of myself. Show me the marker pattern and its relationship to neighboring clusters first; do you want to start with the top genes or the map?

**You**:

What would make you reject or distrust my answer — what would make it useless to you?

**Alex Rios**:

I’d distrust a label based only on **sitting apart on the map**, or on one suggestive gene. It would be useless without the top genes, comparison with neighboring clusters, and some indication that the 163 cells aren’t just a technical or mixed group. Would you give me the marker evidence first, or your proposed identity and how strongly it’s supported?

**You**:

Have you looked at cluster 4's top marker genes yet, or just noticed it sits apart on the UMAP

**Alex Rios**:

Honestly, just the UMAP separation so far—I haven’t examined cluster 4’s top genes yet. So yes, I got ahead of myself; **sitting apart isn’t enough**. Can you compare its top markers with the neighboring clusters and tell me whether they support a known identity or something genuinely distinct?

**You**:

Do you know cluster 4's %mito and gene/count numbers compared to the other clusters

**Alex Rios**:

Yes. Cluster 4 has **1,263 genes/cell** and **2.4% mito**. For comparison: cluster 0: 809/1.8; 1: 850/2.3; 2: 828/2.3; 3: 673/2.1; 5: 1,570/2.0; 6: 350/1.6; 7: 2,363/2.0. Do you want to focus next on the marker genes or the UMAP position?

**You**:

Has anyone considered cluster 4 could be a doublet — two cell types captured together — rather than one new type

**Alex Rios**:

Yes, that’s a real alternative—we haven’t ruled out **a doublet**. Cluster 4 has 1,263 genes/cell and 2.4% mito, but those numbers alone don’t decide it; I’d need to see whether its top genes combine two otherwise separate cell-type programs. Can you check that co-expression pattern?

**You**:

Is there anything about cluster 4 you're taking on faith rather than having actually checked

**Alex Rios**:

Quite a bit. I’ve actually checked only that it has **163 cells**, 1,263 genes/cell, 2.4% mito, and appears separate on the UMAP. I haven’t checked its marker genes, whether those markers match one lineage or two, or whether it resembles a known neighboring population—so the “novel” idea is still just my hypothesis. Which of those should we test first?

**You**:

When and in what format would you like my answer

**Alex Rios**:

By **Friday’s PI handoff**, ideally as a short written call: proposed identity, the key top genes, comparison with neighboring clusters, and whether it looks like one population, a doublet, or an artifact. End with a clear confidence level—don’t just say it’s unusual.

**You**:

The question: Cluster 4 (163 of 2,700 cells) sits apart on the UMAP and hasn't been examined yet. I'm determining whether it's a defensible known cell type or a genuinely novel population — using marker genes, comparison to neighboring clusters, and raw-count co-expression, not UMAP position alone.  The data: pbmc3k.h5ad — 2,700 human PBMCs, 13,714 genes, already processed. ad.X is log-normalized expression (for coloring), ad.layers["counts"] is raw UMI counts, ad.obs["leiden"] has cluster labels 0–7 (unnamed), plus per-cell n_genes, total_counts, pct_mito, and 2-D UMAP coordinates in ad.obsm["X_umap"]. Clusters and the map are already computed.  What I've checked so far on cluster 4: 1,263 genes/cell and 2.4% mito — both mid-range compared to the other clusters (e.g. cluster 6 at 350 genes, cluster 7 at 2,363 genes), so the raw QC numbers alone don't point toward "dying cells" or an obvious doublet. Haven't looked at its top marker genes yet — that's next.  What would satisfy you: A defensible identity (known or novel) backed by top genes, comparison with neighboring clusters, and evidence the 163 cells aren't just a technical/mixed group. You'd distrust a label based only on UMAP position or a single gene.  Deliverable: short written call by Friday — proposed identity, key top genes, neighbor comparison, whether it's one population/doublet/artifact, and a stated confidence level.  Anything I'm missing or got wrong here before I start building?

**Alex Rios**:

The figures and deliverable match what I need. I’d phrase the QC part more neutrally: cluster 4 is **1,263 genes/cell and 2.4% mito**, compared with the other clusters; don’t make the QC conclusion until the marker and raw-count patterns are examined. Start with the top marker genes, or jump straight to the co-expression test?
