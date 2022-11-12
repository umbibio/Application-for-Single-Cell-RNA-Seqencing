import scanpy as sc
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path
import numpy as np

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


# class SelectFromCollection:

#     def __init__(self, ax, compareWindow, c1Label, c2Label, alpha_other=0.3):
#         self.canvas = ax.figure.canvas
#         # self.collection = collection
#         self.alpha_other = alpha_other
#         # self.x = x
#         # self.y = y
#         # self.xys = collection.get_offsets()
#         # self.Npts = len(self.xys)
#         self.compareWindow = compareWindow
#         self.c1Label = c1Label
#         self.c2Label = c2Label

#         # Ensure that we have separate colors for each object
#         # self.fc = collection.get_facecolors()
#         # if len(self.fc) == 0:
#         #     raise ValueError('Collection must have a facecolor')
#         # elif len(self.fc) == 1:
#         #     self.fc = np.tile(self.fc, (self.Npts, 1))

#         self.lasso = LassoSelector(ax, onselect=self.onselect)
#         self.ind = []

#     def onselect(self, verts):
#         path = Path(verts)
#         self.ind = np.nonzero(path.contains_points(self.xys))[0]
#         self.fc[:, -1] = self.alpha_other
#         self.fc[self.ind, -1] = 1
#         # self.collection.set_facecolors(self.fc)
#         self.canvas.draw_idle()
#         for selected in self.xys[self.ind]:
#             selectedArray.append(df.loc[(df[self.x] == selected[0]) & (df[self.y] == selected[1])]['Unnamed: 0'].values[0])
#         # self.compareWindow = ClusterSelectPopup(self.c1Label, self.c2Label)
#         # self.compareWindow.setGeometry(QRect(100, 100, 400, 200))
#         # self.compareWindow.show()

#     def disconnect(self):
#         self.lasso.disconnect_events()
#         self.fc[:, -1] = 1
#         # self.collection.set_facecolors(self.fc)
#         self.canvas.draw_idle()


