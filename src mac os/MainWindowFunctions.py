import scanpy as sc
import pandas as pd

class Cluster:
    ClusterName = ""
    ClusterData = []
    def __str__(self): 
        return self.ClusterName+": "+ str(self.ClusterData)

def createXlsx(adata, grouplen,savefolderpath):
    df = pd.DataFrame()
    for group in range(grouplen):
        if df.empty:
            df = sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')
        else:
            df = pd.concat([df,sc.get.rank_genes_groups_df(adata, group=str(group), key='wilcoxon')])
    df.to_excel(savefolderpath+'/Differential Gene Expression Information.xlsx', index=False)