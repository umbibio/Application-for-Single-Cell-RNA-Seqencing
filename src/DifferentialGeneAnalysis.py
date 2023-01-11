import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as anndata
import ThreadHandling as th
import VisualizationPopup as vp
import MainWindowFunctions as mwf

import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.widgets import LassoSelector
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# from PyQt5 import QtWebEngineWidgets
from PyQt5.QtCore import Qt,pyqtSlot, QRunnable,pyqtSignal,QThreadPool
from PyQt5.QtWidgets import QPushButton, QVBoxLayout, QFormLayout, QHBoxLayout, QLabel, QComboBox,QWidget,QLineEdit,QFrame,QTabWidget,QTabBar,QDesktopWidget,QFileDialog

# import plotly.express as px
# import plotly.graph_objects as go

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
fig = None
plottype = ''

class Window(QWidget):

    def __init__(self, parent, adataFetched):
        super(Window, self).__init__(parent)
        self.parent = parent
        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(lambda index: self.tabs.removeTab(index))
        self.tabs.addTab(self.tab1,"Single cell RNA-Seqencing")
        self.tab1.layout = QHBoxLayout()
        self.tabs.tabBar().setTabButton(0, QTabBar.RightSide,None)

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

        self.Separator = QFrame()
        self.Separator.setFrameShape(QFrame.VLine)
        self.Separator.setLineWidth(1)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.compareWindow = None
        self.displayMessageWindow = None
        # self.browser = QtWebEngineWidgets.QWebEngineView(self)

        self.calcDiffGeneButton = QPushButton("Calculate Differential Gene Expression", self)
        self.calcDiffGeneButton.pressed.connect(self.calcDiffGene)
        self.calcDiffGeneButton.hide()

        self.manualClusteringButton = QPushButton("Manual Clustering", self)
        self.manualClusteringButton.pressed.connect(self.manualClustering)
        self.leidenClusteringButton = QPushButton("Automatic Clustering", self)
        self.leidenClusteringButton.pressed.connect(self.leidenClustering)

        self.c1Label = QLineEdit(Cluster1.ClusterName)
        self.c2Label = QLineEdit(Cluster2.ClusterName)
        self.c1Label.setEnabled(False)
        self.c2Label.setEnabled(False)
        self.c1Label.setAlignment(Qt.AlignCenter)
        self.c2Label.setAlignment(Qt.AlignCenter)
        self.compareButton = QPushButton("Compare", self)
        self.compareButton.pressed.connect(self.compare)
        self.resetButton = QPushButton("Reset", self)
        self.resetButton.pressed.connect(self.reset)
        
        self.outerLayout = QVBoxLayout(self)
        self.attributeLayout = QFormLayout()
        self.graphLayout = QVBoxLayout()
        self.clusteringLayout = QHBoxLayout()
        self.comparisonLayout = QHBoxLayout()
        
        self.attributeLayout.addRow(self.typeLabel,self.typeComboBox)
        self.attributeLayout.addRow(self.gLabel,self.gComboBox)
        self.attributeLayout.addWidget(self.plotButton)

        self.clusteringLayout.addWidget(self.manualClusteringButton)
        self.clusteringLayout.addWidget(self.leidenClusteringButton)

        self.comparisonLayout.addWidget(self.c1Label,1)
        self.comparisonLayout.addWidget(self.c2Label,1)
        self.comparisonLayout.addWidget(self.compareButton,1)
        self.comparisonLayout.addWidget(self.resetButton,1)

        self.tab1.layout.addLayout(self.attributeLayout,1)
        self.tab1.layout.addWidget(self.Separator)
        self.tab1.layout.addLayout(self.graphLayout,2)                                                                                                                                                                                                   
        # layout.setContentsMargins(left, top, right, bottom)
        self.attributeLayout.setContentsMargins(0, 0, 0, 0)
        self.tab1.layout.setContentsMargins(50, 20, 50, 50)
        self.tab1.setLayout(self.tab1.layout)
        self.outerLayout.addWidget(self.tabs)
        self.setLayout(self.outerLayout)
        self.hideClusteringLayout()
        self.hideComparisonLayout()
        self.plotScatter()
        self.threadpool = QThreadPool()

    def reset(self, *args):
        Cluster1.ClusterName = ""
        Cluster1.ClusterData = []
        Cluster2.ClusterName = ""
        Cluster2.ClusterData = []
        self.c1Label.setText("")
        self.c2Label.setText("")

    def hideClusteringLayout(self, *args):
        self.graphLayout.removeItem(self.clusteringLayout)
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
    
    def hideComparisonLayout(self, *args):
        self.graphLayout.removeItem(self.comparisonLayout)
        self.c1Label.hide()
        self.c2Label.hide()
        self.compareButton.hide()
        self.resetButton.hide()
    
    def scatterPlotFigure(self, color):
        global fig
        if plottype == 'PCA':
            fig = sc.pl.pca(adata, color=color, return_fig=True)
        elif plottype == 'UMAP':
            fig = sc.pl.umap(adata, color=color, return_fig=True)
        return fig

    def plotScatter(self, *args):
        global plottype
        self.tabs.removeTab(1)
        self.reset()
        self.hideComparisonLayout()
        self.hideClusteringLayout()
        self.graphLayout.removeWidget(self.calcDiffGeneButton)
        self.calcDiffGeneButton.hide()
        adata.obs['leiden'] =None
        plottype = self.typeComboBox.currentText()
        genetype = self.gComboBox.currentText()
        gene = genes.index[genes['gene_ids'] == genetype].tolist()[0]
        self.figure.clear() 
        self.canvas.close()
        self.scatterPlotFigure(gene)
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        self.graphLayout.addWidget(self.canvas)
        self.graphLayout.addLayout(self.clusteringLayout)
        self.manualClusteringButton.show()
        self.leidenClusteringButton.show()

        # LASSO SELECT USING PLOTLY
        # dfadata = adata.to_df()
        # obsm = adata.obsm.to_df()
        # fig = px.scatter(x = obsm.X_pca1, y = obsm.X_pca2,labels={'x':'PC1', 'y':'PC2'}, title=gene, template="simple_white",color=dfadata[gene])
        # fig.update_layout(paper_bgcolor="#ffffff")
        # f = go.FigureWidget(fig)
        # self.browser.setHtml(f.to_html(include_plotlyjs='cdn'))
        # print(self.browser.page().profile())
        # # scatter = f.data[0]
        # # scatter.on_selection(self.update_point)
        # self.graphLayout.addWidget(self.browser)

    # def update_point(trace, points, selector):
    #     print(points)

    def compare(self, *args):
        global grouplen
        if Cluster1.ClusterData and Cluster2.ClusterData:
            clusterCategory = pd.Series(["0","1"], dtype='category')
            for barcode in Cluster1.ClusterData:
                adata.obs.loc[barcode,'leiden'] = clusterCategory[0]
            for barcode in Cluster2.ClusterData:
                adata.obs.loc[barcode,'leiden'] = clusterCategory[1]
            grouplen = len(adata.obs['leiden'].value_counts())
            self.canvas.close()
            self.scatterPlotFigure('leiden')
            self.canvas = FigureCanvas(fig)
            self.canvas.draw()
            self.hideComparisonLayout()
            self.graphLayout.addWidget(self.canvas,1)
            self.graphLayout.addWidget(self.calcDiffGeneButton)
            self.calcDiffGeneButton.show()

    def leidenClustering(self, *args):
        global AutoClustering, grouplen
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
        sc.tl.leiden(adata)
        grouplen = len(adata.obs['leiden'].value_counts())
        if grouplen>2:
            AutoClustering = True
        self.canvas.close()
        self.scatterPlotFigure('leiden')
        self.canvas = FigureCanvas(fig)
        self.canvas.draw()
        self.graphLayout.addWidget(self.canvas,1)
        self.graphLayout.addWidget(self.calcDiffGeneButton)
        self.calcDiffGeneButton.show()

    def calcDiffGene(self, *args):
        self.calcDiffGeneButton.setEnabled(False)
        thread = calcDiffGeneThread(self)
        thread.signal.progress_signal.connect(self.reportProgress)
        thread.signal.return_signal.connect(self.completePreprocess)
        self.threadpool.start(thread)

    def completePreprocess(self):
        self.calcDiffGeneButton.setEnabled(True)
        self.displayMessageWindow.close()
        self.parent.parent.statusbar.clearMessage()
        self.tab2 = vp.visualizationPopup(grouplen,AutoClustering,adata,self)
        self.tabs.addTab(self.tab2,"Differential Gene Expression Representation")
        self.tabs.setCurrentIndex(1)

    def reportProgress(self,signal):
        print(signal)
        self.parent.parent.statusbar.showMessage(signal)
        self.displayMessageWindow = th.displayMessagePopup(signal)
        self.displayMessageWindow.show()
    
    def manualClustering(self, *args):
        selector = SelectFromCollection(fig.axes[0], fig.axes[0].collections[0], self.c1Label, self.c2Label)
        self.manualClusteringButton.hide()
        self.leidenClusteringButton.hide()
        global obsm, AutoClustering
        AutoClustering=False
        obsm = adata.obsm.to_df()
        self.graphLayout.addLayout(self.comparisonLayout)
        self.c1Label.show()
        self.c2Label.show()
        self.compareButton.show()
        self.resetButton.show()

class calcDiffGeneThread(QRunnable):
    signal = pyqtSignal()

    def __init__(self,parent):
        super(calcDiffGeneThread, self).__init__()     
        self.signal = th.Signals() 
        self.parent=parent   

    @pyqtSlot()
    def run(self):
        self.signal.progress_signal.emit('Executing differential gene analysis...')
        sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='wilcoxon')
        self.signal.progress_signal.emit('Calculating dendrogram...')
        sc.tl.dendrogram(adata, 'leiden')
        self.signal.progress_signal.emit('Creating xlsx file of differential gene expression data...')
        savefolderpath  = QFileDialog.getExistingDirectory(self.parent,'Select folder to save files')+'/CellVisualization'
        if not os.path.isdir(savefolderpath):
            os.mkdir(os.path.normpath(savefolderpath))
        mwf.createXlsx(adata,grouplen,savefolderpath)
        self.signal.return_signal.emit()
        

class SelectFromCollection:

    def __init__(self, ax, collection, c1Label, c2Label, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other
        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)
        self.c1Label = c1Label
        self.c2Label = c2Label
        if plottype == 'PCA':
            self.x='X_pca1'
            self.y='X_pca2'
        elif plottype == 'UMAP':
            self.x='X_umap1'
            self.y='X_umap2'
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
        self.center()
        self.textbox = QLineEdit(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
        self.button = QPushButton('OK', self)
        self.button.setStyleSheet('background-color: #7ad0eb;')
        self.button.move(20,80)
        self.button.resize(280,40)
        self.c1Label = c1Label
        self.c2Label = c2Label
        self.button.clicked.connect(self.on_click)
        grid = QVBoxLayout()
        grid.addWidget(self.textbox)
        grid.addWidget(self.button)

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
    
    def on_click(self):
        global selectedArray
        if selectedArray:
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