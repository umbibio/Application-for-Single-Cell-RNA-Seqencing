import scanpy as sc
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import numpy as np
import pandas as pd

class Cluster:
    ClusterName = ""
    ClusterData = []
    def __str__(self): 
        return self.ClusterName+": "+ str(self.ClusterData)

def preprocessAnnData(adata):
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    adata.write('./PreprocessedData/adata.h5ad')
    return adata


def createXlsx(adata, grouplen):
    df = pd.DataFrame()
    for group in range(grouplen):
        if df.empty:
            df = sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')
        else:
            df.append(sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon'))
    df.to_excel('Differential Gene Expression Information.xlsx', index=False)

#ELEPHANT
#Check which ones can have dendrogram
#NonNAN values check
def generateGraphs(adata,AutoClustering):
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key='wilcoxon', show=False, save='.png')
    sc.pl.dendrogram(adata, 'leiden', show=False, save='.png')
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram = False, show=False, save='.png')
    sc.pl.rank_genes_groups_violin(adata, n_genes=5, key='wilcoxon', show=False, save='.png')
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
    if AutoClustering:
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key='wilcoxon', groupby='leiden', show_gene_labels = True, show=False, save='.png', dendrogram = False)
        sc.pl.rank_genes_groups_tracksplot(adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram= False)
        