import MainWindowFunctions as mwf
import VisualizationPopup as vp
from PyQt5.QtWidgets import QDialog, QPushButton, QVBoxLayout, QFormLayout, QHBoxLayout, QLabel, QComboBox,QWidget,QLineEdit
import pandas as pd
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtCore import QRect,QFile,QTextStream,Qt
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as anndata
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path

AutoClustering = False
Cluster1 = mwf.Cluster()
Cluster2 = mwf.Cluster()
selectedArray = []
gene_ids= []
genes = pd.DataFrame()
obsm = pd.DataFrame()
adata = anndata.AnnData(
  X = None,
  obs = None,
  var = None)
grouplen = 0

class Window(QDialog):

    # constructor
    def __init__(self, adataFetched):
        super(Window, self).__init__()
        self.setWindowTitle('Single cell RNA-Seq')
        self.showMaximized()
        self.setWindowFlag(Qt.WindowMinimizeButtonHint, True)
        self.setWindowFlag(Qt.WindowMaximizeButtonHint, True)

        global adata, genes, gene_ids
        adata = adataFetched
        genes = pd.DataFrame(adata.var.gene_ids)
        gene_ids = genes['gene_ids'].to_numpy()
        self.typeComboBox = QComboBox(self)
        self.typeComboBox.addItems(['PCA', 'UMAP'])
        self.typeLabel = QLabel("Type:")
        self.typeLabel.setBuddy(self.typeComboBox)
        
        self.gComboBox = QComboBox()
        self.gComboBox.addItems(gene_ids)
        self.gLabel = QLabel("Gene:")
        self.gLabel.setBuddy(self.gComboBox)

        self.plotButton = QPushButton("Plot Current Attributes", self)
        self.plotButton.pressed.connect(self.plotScatter)

        self.manualClusteringButton = QPushButton("Manual Clustering", self)
        self.manualClusteringButton.pressed.connect(self.manualClustering)
        self.leidenClusteringButton = QPushButton("Automatic Clustering", self)
        self.leidenClusteringButton.pressed.connect(self.leidenClustering)
        self.calcDiffGeneButton = QPushButton("Calculate Differential Gene Expression", self)
        self.calcDiffGeneButton.pressed.connect(self.calcDiffGene)
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
        self.calcDiffGeneButton.hide()
        
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.compareWindow = None

        self.c1Label = QLineEdit(Cluster1.ClusterName)
        self.c2Label = QLineEdit(Cluster2.ClusterName)
        self.c1Label.setEnabled(False)
        self.c2Label.setEnabled(False)
        self.c1Label.setAlignment(Qt.AlignCenter)
        self.c2Label.setAlignment(Qt.AlignCenter)
        self.compareButton = QPushButton("Compare", self)
        self.compareButton.pressed.connect(self.compare)
        self.compareButton.hide()

        # self.resetButton = QPushButton("Reset", self)
        # self.resetButton.pressed.connect(self.reset):

        self.outerLayout = QVBoxLayout()
        self.attributeLayout = QFormLayout()
        self.clusteringLayout = QHBoxLayout()
        self.comparisonLayout = QHBoxLayout()
        

        self.attributeLayout.addRow(self.typeLabel,self.typeComboBox)
        self.attributeLayout.addRow(self.gLabel,self.gComboBox)

        self.clusteringLayout.addWidget(self.manualClusteringButton)
        self.clusteringLayout.addWidget(self.leidenClusteringButton)

        self.comparisonLayout.addWidget(self.c1Label,1)
        self.comparisonLayout.addWidget(self.c2Label,1)
        self.comparisonLayout.addWidget(self.compareButton,1)

        self.outerLayout.addLayout(self.attributeLayout)
        self.outerLayout.addWidget(self.plotButton)
        self.outerLayout.addWidget(self.canvas,1)
        
        # self.outerLayout.addLayout(self.comparisonLayout)
        

        # layout.setContentsMargins(left, top, right, bottom)
        self.attributeLayout.setContentsMargins(0, 0, 0, 0)
        self.outerLayout.setContentsMargins(50, 20, 50, 50)
        self.setLayout(self.outerLayout)

    # def reset(self, *args):
    #     Cluster1.ClusterName = ""
    #     Cluster1.ClusterData = []
    #     Cluster2.ClusterName = ""
    #     Cluster2.ClusterData = []
    #     self.c1Label.setText("")
    #     self.c2Label.setText("")

    def plotScatter(self, *args):
        
        self.outerLayout.removeItem(self.clusteringLayout)
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
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
            
        selector = SelectFromCollection(fig.axes[0], fig.axes[0].collections[0], plottype, self.c1Label, self.c2Label)
        
        self.canvas.draw()
        self.outerLayout.addWidget(self.canvas,1)
        self.outerLayout.addLayout(self.clusteringLayout)
        self.manualClusteringButton.show()
        self.leidenClusteringButton.show()
        # self.outerLayout.addWidget(self.resetButton)
        

    def compare(self, *args):
        # self.outerLayout.removeWidget(self.resetButton)
        global grouplen
        clusterCategory = pd.Series(["0","1"], dtype='category')
        plottype = self.typeComboBox.currentText()
        for barcode in Cluster1.ClusterData:
            adata.obs.loc[barcode,'leiden'] = clusterCategory[0]
        for barcode in Cluster2.ClusterData:
            adata.obs.loc[barcode,'leiden'] = clusterCategory[1]
        grouplen = len(adata.obs['leiden'].value_counts())
        self.canvas.close()
        if plottype == 'PCA':
            fig = sc.pl.pca(adata, color='leiden', return_fig=True)
            # self.figure.add_axes(ax)
        elif plottype == 'UMAP':
            fig = sc.pl.umap(adata, color='leiden', return_fig=True)
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        self.outerLayout.removeItem(self.comparisonLayout)
        self.c1Label.hide()
        self.c2Label.hide()
        self.compareButton.hide()
        self.outerLayout.addWidget(self.canvas,1)
        self.outerLayout.addWidget(self.calcDiffGeneButton)
        self.calcDiffGeneButton.show()
        # self.outerLayout.addWidget(self.resetButton)


    def leidenClustering(self, *args):
        global AutoClustering, grouplen
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
        plottype = self.typeComboBox.currentText()
        sc.tl.leiden(adata)
        grouplen = len(adata.obs['leiden'].value_counts())
        if grouplen>2:
            AutoClustering = True
        self.canvas.close()
        if plottype == 'PCA':
            fig = sc.pl.pca(adata, color='leiden', return_fig=True)
            # self.figure.add_axes(ax)
        elif plottype == 'UMAP':
            fig = sc.pl.umap(adata, color='leiden', return_fig=True)
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        self.outerLayout.addWidget(self.canvas,1)
        self.outerLayout.addWidget(self.calcDiffGeneButton)
        self.calcDiffGeneButton.show()
        # self.outerLayout.addWidget(self.resetButton)

    def calcDiffGene(self, *args):
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
        sc.tl.dendrogram(adata, 'leiden')
        mwf.createXlsx(adata,grouplen)
        mwf.generateGraphs(adata, AutoClustering)
        self.visualizationWindow = vp.visualizationPopup(grouplen,AutoClustering)
        self.visualizationWindow.show()
        self.canvas.close()
        self.outerLayout.removeWidget(self.calcDiffGeneButton)
        self.calcDiffGeneButton.hide()
        plottype = self.typeComboBox.currentText()
        if plottype == 'PCA':
            fig = sc.pl.pca(adata, color='leiden', return_fig=True)
            # self.figure.add_axes(ax)
        elif plottype == 'UMAP':
            fig = sc.pl.umap(adata, color='leiden', return_fig=True)
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        self.outerLayout.addWidget(self.canvas,1)
        self.outerLayout.addWidget(self.calcDiffGeneButton)
        self.calcDiffGeneButton.show()
    
    def manualClustering(self, *args):
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
        # self.outerLayout.removeWidget(self.resetButton)
        global obsm
        obsm = adata.obsm.to_df()
        # self.clusteringLayout.setVisible(False)
        self.outerLayout.addLayout(self.comparisonLayout)
        self.compareButton.show()
        # self.outerLayout.addWidget(self.resetButton)


class SelectFromCollection:

    def __init__(self, ax, collection, plottype, c1Label, c2Label, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other
        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)
        self.plottype = plottype
        self.c1Label = c1Label
        self.c2Label = c2Label

        if plottype == 'PCA':
            self.x='X_pca1'
            self.y='X_pca2'
        elif plottype == 'UMAP':
            self.x='X_umap1'
            self.y='X_umap2'

        # Ensure that we have separate colors for each object
        self.fc = collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Npts, 1))
        self.lasso = LassoSelector(ax, onselect=self.onselect)
        self.ind = []

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()
        for selected in self.xys[self.ind]:
            selectedGene = obsm.index[(obsm[self.x] == selected[0]) & (obsm[self.y] == selected[1])][0]
            selectedArray.append(selectedGene)
        self.compareWindow = ClusterSelectPopup(self.c1Label, self.c2Label)
        # self.compareWindow.setGeometry(QRect(100, 100, 400, 200))
        self.compareWindow.show()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()


class ClusterSelectPopup(QWidget):
    def __init__(self, c1Label, c2Label):
        QWidget.__init__(self)
        self.setStyleSheet('background-color:white;')
        self.setWindowTitle('Cluster Name')
        self.textbox = QLineEdit(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
        # self.textbox.setAlignment(Qt.AlignCenter)
        self.button = QPushButton('OK', self)
        self.button.setStyleSheet('background-color: #7ad0eb;')
        self.button.move(20,80)
        self.button.resize(280,40)
        self.c1Label = c1Label
        self.c2Label = c2Label
        # connect button to function on_click
        self.button.clicked.connect(self.on_click)
        grid = QVBoxLayout()
        grid.addWidget(self.textbox)
        grid.addWidget(self.button)

    def on_click(self):
        global selectedArray
        if not Cluster1.ClusterName:
            Cluster1.ClusterName = self.textbox.text()
            Cluster1.ClusterData = selectedArray
            self.c1Label.setText(Cluster1.ClusterName)
        else :
            Cluster2.ClusterName = self.textbox.text()
            Cluster2.ClusterData = selectedArray
            self.c2Label.setText(Cluster2.ClusterName)
        selectedArray=[]
        self.close()