import MainWindowFunctions as mwf
import VisualizationPopup as vp
from PyQt5.QtWidgets import QDialog, QPushButton, QVBoxLayout, QFormLayout, QHBoxLayout, QLabel, QComboBox
from PyQt5.QtCore import QRect
import pandas as pd
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as anndata

AutoClustering = False
Cluster1 = mwf.Cluster()
Cluster2 = mwf.Cluster()
gene_ids= []
genes = pd.DataFrame()
adata = anndata.AnnData(
  X = None,
  obs = None,
  var = None)

class Window(QDialog):

    # constructor
    def __init__(self, adataFetched):
        super(Window, self).__init__()
        self.setWindowTitle('Single cell RNA-Seq')
        self.showMaximized()

        global adata, genes, gene_ids
        adata = adataFetched
        genes = pd.DataFrame(adata.var.gene_ids)
        gene_ids = genes['gene_ids'].to_numpy()
        self.typeComboBox = QComboBox(self)
        self.typeComboBox.addItems(['PCA', 'UMAP'])
        self.typeLabel = QLabel("Type:")
        self.typeLabel.setBuddy(self.typeComboBox)
        self.typeComboBox.resize(165, self.typeComboBox.height())
        
        self.gComboBox = QComboBox()
        self.gComboBox.addItems(gene_ids)
        self.gLabel = QLabel("Gene:")
        self.gLabel.setBuddy(self.gComboBox)

        self.plotButton = QPushButton("Plot Current Attributes", self)
        self.plotButton.pressed.connect(self.plotScatter)

        self.leidenClusteringButton = QPushButton("Automatic Clustering", self)
        self.leidenClusteringButton.pressed.connect(self.leidenClustering)
        self.calcDiffGeneButton = QPushButton("Calculate Differential Gene Expression", self)
        self.calcDiffGeneButton.pressed.connect(self.calcDiffGene)
        

        self.plotComboBox = QComboBox()
        self.plotComboBox.addItems(['Histogram','Differential gene expression','Dot Plot', 'Violin', 'Stacked Violin', 'Matrix Plot', 'Heatmap', 'Tracksplot'])
        self.plotLabel = QLabel("Plot Type:")
        self.plotLabel.setBuddy(self.plotComboBox)
        # self.plotComboBox.currentIndexChanged.connect(self.updateGraph)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.compareWindow = None

        self.c1Label = QLabel(Cluster1.ClusterName)
        self.c2Label = QLabel(Cluster2.ClusterName)
        self.c1Label.setStyleSheet("border: 1px solid black;")
        self.c2Label.setStyleSheet("border: 1px solid black;")
        self.compareButton = QPushButton("Compare", self)
        # self.compareButton.pressed.connect(self.compare)

        self.resetButton = QPushButton("Reset", self)
        # self.resetButton.pressed.connect(self.reset)

        self.outerLayout = QVBoxLayout()
        self.attributeLayout = QFormLayout()
        self.plotSelectionLayout = QFormLayout()
        self.comparisonLayout = QHBoxLayout()

        self.attributeLayout.addRow(self.typeLabel,self.typeComboBox)
        self.attributeLayout.addRow(self.gLabel,self.gComboBox)

        self.plotSelectionLayout.addRow(self.plotLabel,self.plotComboBox)

        self.comparisonLayout.addWidget(self.c1Label)
        self.comparisonLayout.addWidget(self.c2Label)
        self.comparisonLayout.addWidget(self.compareButton)

        self.outerLayout.addLayout(self.attributeLayout)
        self.outerLayout.addWidget(self.plotButton )
        
        # self.outerLayout.addLayout(self.plotSelectionLayout)
        # self.outerLayout.addWidget(self.canvas)
        # self.outerLayout.addLayout(self.comparisonLayout)
        

        # layout.setContentsMargins(left, top, right, bottom)
        self.attributeLayout.setContentsMargins(0, 0, 0, 0)
        self.outerLayout.setContentsMargins(50, 20, 50, 50)
        self.setLayout(self.outerLayout)

    def reset(self, *args):
        Cluster1.ClusterName = ""
        Cluster1.ClusterData = []
        Cluster2.ClusterName = ""
        Cluster2.ClusterData = []
        self.c1Label.setText("")
        self.c2Label.setText("")

    def plotScatter(self, *args):
        
        # finding the content of current item in combo box
        plottype = self.typeComboBox.currentText()
        genetype = self.gComboBox.currentText()
        gene = genes.index[genes['gene_ids'] == genetype].tolist()[0]
        self.figure.clear() #clear current figure
        self.canvas.close()
        if plottype == 'PCA':
            fig = sc.pl.pca(adata, color=gene, return_fig=True)
            self.canvas = FigureCanvas(fig)
            # self.figure.add_axes(ax)
        elif plottype == 'UMAP':
            fig = sc.pl.umap(adata, color=gene, return_fig=True)
            self.canvas = FigureCanvas(fig)
            
        # selector = mwf.SelectFromCollection(fig.axes, fig.axes.collections[0], self.compareWindow, self.c1Label, self.c2Label)
        
        self.canvas.draw()
        self.outerLayout.addWidget(self.canvas)
        self.outerLayout.addWidget(self.leidenClusteringButton)
        self.outerLayout.addWidget(self.resetButton)
        # draw new graph
        # self.canvas.draw()

    def leidenClustering(self, *args):
        AutoClustering = True
        self.outerLayout.removeWidget(self.leidenClusteringButton)
        self.outerLayout.removeWidget(self.resetButton)
        plottype = self.typeComboBox.currentText()
        sc.tl.leiden(adata)
        self.canvas.close()
        if plottype == 'PCA':
            fig = sc.pl.pca(adata, color='leiden', return_fig=True)
            # self.figure.add_axes(ax)
        elif plottype == 'UMAP':
            fig = sc.pl.umap(adata, color='leiden', return_fig=True)
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        self.outerLayout.addWidget(self.canvas)
        self.outerLayout.addWidget(self.calcDiffGeneButton)
        self.outerLayout.addWidget(self.resetButton)

    def calcDiffGene(self, *args):
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
        mwf.createXlsx(adata)
        sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, key='wilcoxon', show=False, save='.png')
        sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key='wilcoxon', groupby='leiden', dendrogram = False, show=False, save='.png')
        sc.pl.rank_genes_groups_violin(adata, n_genes=5, key='wilcoxon', show=False, save='.png')
        sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
        sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
        sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, key='wilcoxon', groupby='leiden', show_gene_labels = True, show=False, save='.png', dendrogram = False)
        sc.pl.rank_genes_groups_tracksplot(adata, n_genes=5, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram= False)
        self.visualizationWindow = vp.visualizationPopup()
        # self.visualizationWindow.setGeometry(QRect(100, 100, 1000, 800))
        self.visualizationWindow.show()
