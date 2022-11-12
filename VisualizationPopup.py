import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.backends.backend_qt5agg import \
    FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QComboBox, QLabel, QVBoxLayout, QWidget


class visualizationPopup(QWidget):
    def __init__(self, adata, AutoClustering):
        QWidget.__init__(self)
        self.setWindowTitle('Cell Visualization App')
        self.adata = adata
        self.AutoClustering = AutoClustering
        self.plotComboBox = QComboBox()
        self.plotComboBox.addItems(['','Dot Plot', 'Violin', 'Stacked Violin', 'Matrix Plot', 'Heatmap', 'Tracksplot'])
        self.plotLabel = QLabel("Plot Type:")
        self.plotLabel.setBuddy(self.plotComboBox)
        self.plotComboBox.currentIndexChanged.connect(self.updateGraph())

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        grid = QVBoxLayout()
        grid.addWidget(self.plotLabel)
        grid.addWidget(self.plotComboBox)
        grid.addWidget(self.canvas)

        self.setLayout(grid)
        

    def updateGraph(self, *args):

        self.figure.clear() #clear current figure
        self.canvas.close()
        plot = self.plotComboBox.currentText();
        # if plot == '':
        #     sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key='wilcoxon')
        # el
        if  plot == 'Dot Plot':
            fig = sc.pl.rank_genes_groups_dotplot(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram = self.AutoClustering, return_fig=True)
        elif plot == 'Violin':
            fig = sc.pl.rank_genes_groups_violin(self.adata, groups='1', n_genes=5, key='wilcoxon', return_fig=True)
        elif plot == 'Stacked Violin':
            fig = sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram = self.AutoClustering, return_fig=True)
        elif plot == 'Matrix Plot':
            fig = sc.pl.rank_genes_groups_matrixplot(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram = self.AutoClustering, return_fig=True)
        elif plot == 'Heatmap':
            fig = sc.pl.rank_genes_groups_heatmap(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', show_gene_labels = True, dendrogram = self.AutoClustering, return_fig=True)
        elif plot == 'Tracksplot':
            # self.figure.clear() #clear current figure
            # ax = sc.pl.tracksplot(adata, marker_genes_dict, groupby='leiden', dendrogram=True)
            fig = sc.pl.rank_genes_groups_tracksplot(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram= self.AutoClustering, return_fig=True)
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()