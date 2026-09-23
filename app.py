from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scanpy as sc
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse

DATA_PATH = Path(__file__).parent / "data" / "pbmc3k.h5ad"
adata = None
markers: dict[str, list[dict[str, Any]]] = {}


def scalar(value):
    return value.item() if hasattr(value, "item") else value


def values_for_gene(gene: str) -> np.ndarray:
    if gene not in adata.var_names:
        raise HTTPException(404, f"Unknown gene: {gene}")
    col = adata[:, gene].X
    return np.asarray(col.toarray()).ravel() if hasattr(col, "toarray") else np.asarray(col).ravel()


@asynccontextmanager
async def lifespan(app: FastAPI):
    global adata, markers
    adata = sc.read_h5ad(DATA_PATH)
    # Rank markers once; the dataset's existing leiden labels are used unchanged.
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", n_genes=20)
    groups = list(adata.obs["leiden"].cat.categories if hasattr(adata.obs["leiden"], "cat") else sorted(adata.obs["leiden"].unique()))
    for group in groups:
        markers[str(group)] = []
        for i in range(min(20, len(adata.uns["rank_genes_groups"]["names"]))):
            name = scalar(adata.uns["rank_genes_groups"]["names"][i][str(group)])
            score = scalar(adata.uns["rank_genes_groups"]["scores"][i][str(group)])
            logfc = scalar(adata.uns["rank_genes_groups"].get("logfoldchanges", np.full_like(adata.uns["rank_genes_groups"]["scores"], np.nan))[i][str(group)])
            pval = scalar(adata.uns["rank_genes_groups"].get("pvals_adj", np.full_like(adata.uns["rank_genes_groups"]["scores"], np.nan))[i][str(group)])
            markers[str(group)].append({"gene": name, "score": float(score), "logfoldchange": float(logfc), "pval_adj": float(pval)})
    yield


app = FastAPI(title="PBMC3k Explorer", lifespan=lifespan)


@app.get("/", response_class=HTMLResponse)
async def index():
    return HTML


@app.get("/api/umap")
async def umap():
    coords = np.asarray(adata.obsm["X_umap"])
    return {"x": coords[:, 0].tolist(), "y": coords[:, 1].tolist(), "cluster": adata.obs["leiden"].astype(str).tolist()}


@app.get("/api/genes")
async def genes():
    return {"genes": adata.var_names.astype(str).tolist()}


@app.get("/api/expression/{gene}")
async def expression(gene: str):
    return {"gene": gene, "values": values_for_gene(gene).astype(float).tolist()}


@app.get("/api/quality")
async def quality():
    fields = [f for f in ("n_genes", "total_counts", "pct_mito") if f in adata.obs]
    return {field: pd.to_numeric(adata.obs[field]).tolist() for field in fields}


@app.get("/api/clusters")
async def clusters():
    result = {}
    for cluster in sorted(adata.obs["leiden"].astype(str).unique(), key=lambda x: int(x)):
        subset = adata.obs[adata.obs["leiden"].astype(str) == cluster]
        result[cluster] = {
            "n_cells": len(subset),
            "quality": {field: float(pd.to_numeric(subset[field]).mean()) for field in ("n_genes", "total_counts", "pct_mito") if field in subset},
            "markers": markers[cluster],
        }
    return result


@app.get("/api/analysis")
async def analysis():
    coords = np.asarray(adata.obsm["X_umap"])
    labels = adata.obs["leiden"].astype(str).to_numpy()
    clusters_sorted = sorted(np.unique(labels), key=lambda x: int(x))
    centroids = {c: coords[labels == c].mean(axis=0) for c in clusters_sorted}
    distances = {c: {d: float(np.linalg.norm(centroids[c] - centroids[d])) for d in clusters_sorted if d != c} for c in clusters_sorted}
    nearest_4 = [{"cluster": c, "distance": d} for c, d in sorted(distances["4"].items(), key=lambda item: item[1])]
    return {"centroids": {c: centroids[c].tolist() for c in clusters_sorted}, "distances": distances, "cluster_4_nearest": nearest_4, "clusters": markers}


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PBMC3k Explorer</title><script src="https://cdn.tailwindcss.com"></script><script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script></head>
<body class="bg-slate-50 text-slate-900"><main class="mx-auto max-w-7xl p-6 sm:p-8 lg:p-10"><header class="mb-8"><h1 class="text-3xl font-bold sm:text-4xl">PBMC3k Explorer</h1><p class="mt-2 text-base text-slate-600">2,700 cells · 13,714 genes · explore clusters, expression and quality</p></header>
<section class="grid gap-8 lg:grid-cols-[19rem_1fr]"><aside class="space-y-6 rounded-xl bg-white p-6 shadow-sm"><label class="block text-base font-medium">Colour by<select id="mode" class="mt-2 w-full rounded border p-3 text-base"><option value="cluster">Leiden cluster</option><option value="gene">Gene expression</option><option value="quality">Quality metric</option></select></label><label id="geneWrap" class="hidden text-base font-medium">Search gene<input id="geneSearch" type="search" autocomplete="off" placeholder="Type a gene symbol…" class="mt-2 w-full rounded border p-3 text-base"><select id="gene" size="6" class="mt-2 w-full rounded border p-3 text-base"></select></label><label id="qualityWrap" class="hidden text-base font-medium">Metric<select id="quality" class="mt-2 w-full rounded border p-3 text-base"><option>n_genes</option><option>total_counts</option><option>pct_mito</option></select></label><div id="stats" class="text-base text-slate-600"></div></aside><div class="min-w-0 rounded-xl bg-white p-5 shadow-sm sm:p-7"><div id="plot" class="h-[60vh] min-h-[420px] w-full"></div></div></section>
<section class="mt-8 rounded-xl bg-white p-6 shadow-sm sm:p-8"><h2 class="text-2xl font-semibold">Cluster top marker genes</h2><p class="mt-2 text-base text-slate-600">Top markers ranked by score; log-fold-change and adjusted p-value are shown for comparison.</p><div id="clusters" class="mt-6 overflow-x-auto"></div></section>
<section class="mt-10 space-y-8"><div class="rounded-xl bg-white p-6 shadow-sm sm:p-8"><h2 class="text-3xl font-bold">Cluster Identification &amp; Analysis</h2><p class="mt-2 text-base text-slate-600">Computed UMAP centroid distances for cluster 4, ranked from nearest to farthest.</p><ol id="neighbors" class="mt-5 grid gap-3 sm:grid-cols-2 lg:grid-cols-3"></ol></div><div id="analysis" class="space-y-6"></div><p class="px-2 pb-4 text-sm text-slate-500">Gene-to-cell-type interpretation in this analysis was done with assistance from Claude (Anthropic); all marker gene rankings, statistics, and QC numbers above are computed directly from the dataset by this app.</p></section></main>
<script>
const state={}; const $=id=>document.getElementById(id);
async function load(){[state.umap,state.genes,state.quality,state.clusters,state.analysis]=await Promise.all(['/api/umap','/api/genes','/api/quality','/api/clusters','/api/analysis'].map(x=>fetch(x).then(r=>r.json()))); filterGenes(); renderClusters(); renderAnalysis(); update();}
function update(){let mode=$('mode').value; $('geneWrap').classList.toggle('hidden',mode!=='gene');$('qualityWrap').classList.toggle('hidden',mode!=='quality'); let color, title;
if(mode==='cluster'){color=state.umap.cluster;title='Leiden cluster'} else if(mode==='quality'){color=state.quality[$('quality').value];title=$('quality').value}else{let gene=$('gene').value;fetch('/api/expression/'+encodeURIComponent(gene)).then(r=>r.json()).then(x=>draw(x.values,gene));return} draw(color,title)}
function draw(color,title){Plotly.react('plot',[{x:state.umap.x,y:state.umap.y,mode:'markers',type:'scattergl',text:state.umap.cluster,marker:{color:color,colorscale:'Viridis',size:6,opacity:.8,colorbar:{title:title}},hovertemplate:'cluster %{text}<br>x %{x:.2f}, y %{y:.2f}<extra></extra>'}],{margin:{l:45,r:20,t:15,b:45},xaxis:{title:'UMAP 1'},yaxis:{title:'UMAP 2'},paper_bgcolor:'white',plot_bgcolor:'white'})}
function renderClusters(){ $('clusters').innerHTML=`<table class="w-full min-w-[1100px] border-collapse text-base"><thead><tr class="border-b-2 text-left"><th class="p-4">Cluster</th><th class="p-4">Cells</th><th class="p-4">Top marker genes (ranked by score)</th><th class="p-4">Mean QC</th></tr></thead><tbody>${Object.entries(state.clusters).map(([c,v])=>`<tr class="border-b align-top hover:bg-slate-50"><td class="p-4 text-lg font-semibold">${c}</td><td class="p-4">${v.n_cells}</td><td class="p-4"><div class="grid grid-cols-2 gap-x-8 gap-y-2 lg:grid-cols-4">${v.markers.slice(0,8).map((x,i)=>`<div><span class="font-mono font-semibold">${i+1}. ${x.gene}</span><div class="text-sm text-slate-600">LFC ${x.logfoldchange.toFixed(2)} · score ${x.score.toFixed(2)} · p ${x.pval_adj.toExponential(1)}</div></div>`).join('')}</div></td><td class="whitespace-nowrap p-4 text-sm text-slate-600">genes ${v.quality.n_genes?.toFixed(0)}<br>counts ${v.quality.total_counts?.toFixed(0)}<br>mito ${v.quality.pct_mito?.toFixed(1)}%</td></tr>`).join('')}</tbody></table>`}
function renderAnalysis(){const descriptions={0:['CD4 T cells','CD3D and ribosomal markers support T-cell lineage and a resting state.'],1:['CD14+ classical monocytes','LYZ, S100A8/S100A9, FCN1 and TYROBP define classical monocytes.'],2:['NK / cytotoxic T cells','NKG7, GZMA, CST7, CTSW and CCL5 form a cytotoxic effector signature.'],3:['B cells','CD79A/B, MS4A1 and HLA-II genes support B-cell identity.'],4:['FCGR3A+ (CD16+) non-classical monocytes','FCGR3A, FCER1G, LST1 and AIF1 form a coherent myeloid program without T-, B- or NK-cell markers. This is a known, distinct monocyte population, not a novel type.'],5:['Dendritic cells','Concentrated HLA-II, CD74 and FCER1A expression supports dendritic-cell identity.'],6:['Megakaryocytes / platelets','PF4 and PPBP are unmistakable platelet markers in this tiny population.'],7:['Likely artifact / doublet — no cell-type assignment','Housekeeping and cell-cycle genes lack a lineage-defining program; weak statistics and unusually high RNA content in few cells flag a likely doublet or technical artifact.']}; $('neighbors').innerHTML=state.analysis.cluster_4_nearest.map((n,i)=>`<li class="rounded-lg border p-4 text-base"><b>${i+1}. Cluster ${n.cluster}</b><span class="ml-2 text-slate-600">distance ${n.distance.toFixed(3)}</span></li>`).join(''); $('analysis').innerHTML=Object.entries(state.clusters).map(([c,v])=>{let [identity,reason]=descriptions[c];let evidence=v.markers.slice(0,8).map(x=>`${x.gene} (LFC ${x.logfoldchange.toFixed(2)}, score ${x.score.toFixed(2)}, p ${x.pval_adj.toExponential(1)})`).join('; ');let qc=`${v.quality.n_genes?.toFixed(0)} genes, ${v.quality.total_counts?.toFixed(0)} counts, ${v.quality.pct_mito?.toFixed(1)}% mito`;return `<article class="rounded-xl ${c==='4'?'border-4 border-indigo-400 bg-indigo-50 p-8':'bg-white p-6'} shadow-sm"><h3 class="text-xl font-bold">Cluster ${c}: ${identity}</h3><p class="mt-4 text-base leading-8">${reason}</p><p class="mt-3 text-sm leading-7 text-slate-700"><b>Computed evidence:</b> ${evidence}. <b>QC:</b> ${qc}.</p>${c==='4'?'<p class="mt-4 text-base font-semibold">Verdict: high-confidence FCGR3A+ non-classical monocytes; closest centroid neighbor is shown above.</p>':''}</article>`}).join('')}
function filterGenes(){const query=$('geneSearch').value.trim().toLowerCase();$('gene').innerHTML='';state.genes.genes.filter(g=>g.toLowerCase().includes(query)).slice(0,100).forEach(g=>$('gene').add(new Option(g,g)));if(!$('gene').value&&$('gene').options.length)$('gene').selectedIndex=0}
$('mode').onchange=update;$('gene').onchange=update;$('geneSearch').oninput=()=>{filterGenes();update()};$('quality').onchange=update;load();
</script></body></html>'''
