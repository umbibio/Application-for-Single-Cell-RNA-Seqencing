import cv2 as cv
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.backends.backend_qt5agg import \
    FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QComboBox, QLabel, QVBoxLayout, QWidget


class visualizationPopup(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.setWindowTitle('Differential Gene Expression Representation')
        self.showMaximized()
        self.plotComboBox = QComboBox()
        self.plotComboBox.addItems(['Gene Ranking','Dot Plot', 'Violin', 'Stacked Violin', 'Matrix Plot', 'Heatmap', 'Tracksplot'])
        self.plotLabel = QLabel("Plot Type:")
        self.plotLabel.setBuddy(self.plotComboBox)
        self.plotComboBox.currentIndexChanged.connect(self.updateGraph)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        grid = QVBoxLayout()
        grid.addWidget(self.plotLabel)
        grid.addWidget(self.plotComboBox)
        grid.addWidget(self.canvas)

        self.setLayout(grid)
        

    def updateGraph(self, *args):

        # self.canvas.close()
        plot = self.plotComboBox.currentText()
        # if plot == '':
        #     sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key='wilcoxon')
        # el
        self.figure.clear() #clear current figure
        self.canvas.flush_events()
        ax = self.figure.add_subplot(111) # create an axis
        if plot == 'Gene Ranking':
            # sc.pl.rank_genes_groups(self.adata, n_genes=25, sharey=False, key='wilcoxon', show=False, save='.png')
            ax.imshow('./figures/dotplot_.png')
        elif  plot == 'Dot Plot':
            # sc.pl.rank_genes_groups_dotplot(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram = False, show=False, save='.png')
            ax.imshow('./figures/dotplot_.png')
        # elif plot == 'Violin':
        #     print()
        #     # ax = self.figure.add_subplot(121) # create an axis
        #     # sc.pl.rank_genes_groups_violin(self.adata, groups='1', n_genes=5, key='wilcoxon', show=False, save='.png')
        #     # ax1 = self.figure.add_subplot(122) # create an axis
        #     # sc.pl.rank_genes_groups_violin(self.adata, groups='2', n_genes=5, key='wilcoxon', ax = ax1)
        # elif plot == 'Stacked Violin':
        #     sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
        # elif plot == 'Matrix Plot':
        #     sc.pl.rank_genes_groups_matrixplot(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
        # elif plot == 'Heatmap':
        #     sc.pl.rank_genes_groups_heatmap(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', show_gene_labels = True, show=False, save='.png', dendrogram = False)
        # elif plot == 'Tracksplot':
        #     # self.figure.clear() #clear current figure
        #     # ax = sc.pl.tracksplot(adata, marker_genes_dict, groupby='clusters', dendrogram=True)
        #     sc.pl.rank_genes_groups_tracksplot(self.adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram= False)
        self.canvas.draw()