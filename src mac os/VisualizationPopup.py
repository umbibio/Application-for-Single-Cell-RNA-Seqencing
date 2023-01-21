import os
import cv2 as cv
import scanpy as sc
import ThreadHandling as th
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from PyQt5.QtCore import Qt,pyqtSlot, QRunnable,pyqtSignal,QThreadPool
from PyQt5.QtWidgets import QComboBox, QLabel, QVBoxLayout, QWidget, QFormLayout,QLineEdit,QPushButton,QDialog,QSlider

n_genes=5
plot=''
basedir = os.path.dirname(__file__)
sc.settings.figdir=basedir+'/figures/'

class visualizationPopup(QWidget):
    def __init__(self,grouplen,AutoClustering,adata,parent):
        QWidget.__init__(self)
        self.threadpool = QThreadPool()
        self.parent=parent
        self.adata = adata
        self.grouplen=grouplen
        self.AutoClustering = AutoClustering
        self.setWindowTitle('Differential Gene Expression Representation')
        self.showMaximized()
        self.plotComboBox = QComboBox()
        self.plotComboBox.addItems(['Dendrogram','Gene Ranking', 'Dot Plot', 'Violin', 'Stacked Violin', 'Matrix Plot', 'Heatmap', 'Tracksplot'])

        self.plotLabel = QLabel("Plot Type:")
        self.plotLabel.setBuddy(self.plotComboBox)
        self.plotComboBox.currentIndexChanged.connect(self.updateGraph)
        self.figure = plt.figure()
        self.canvasPop = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvasPop, self)
        self.plotSelectionLayout = QFormLayout()
        self.plotSelectionLayout.addRow(self.plotLabel,self.plotComboBox)
        self.NGenesWindow = None
        self.displayMessageWindow = None

        self.slider = QSlider(Qt.Horizontal,self)
        # self.slider.setGeometry(60,60,self.grouplen,20)
        self.slider.valueChanged[int].connect(self.updateViolin)
        self.slider.setMaximum(self.grouplen-1)

        self.grid = QVBoxLayout()
        self.grid.addLayout(self.plotSelectionLayout)
        self.grid.addWidget(self.toolbar)
        self.grid.addWidget(self.canvasPop,1)
        self.setLayout(self.grid)
        self.updateGraph()
    
    def updateViolin(self, value):
        self.figure.clear()
        self.canvasPop.flush_events()
        ax = self.figure.add_subplot(111) 
        img = cv.imread(os.path.join(basedir,'figures','rank_genes_groups_leiden_'+str(value)+'.png'))
        ax.imshow(img)
        self.canvasPop.draw()

        
    def updateGraph(self, *args):
        global plot
        self.grid.removeWidget(self.slider)
        plot = self.plotComboBox.currentText()
        self.figure.clear()
        self.canvasPop.flush_events()
        if plot == 'Dendrogram':
            self.parent.parent.parent.statusbar.showMessage("Generating Dendrogram...")
            sc.pl.dendrogram(self.adata, 'leiden', show=False, save='.png')
            ax = self.figure.add_subplot(111)
            img = cv.imread(os.path.join(basedir,'figures','dendrogram.png'))
            ax.imshow(img)
            self.canvasPop.draw()
        else:
            self.plotComboBox.setEnabled(False)
            self.NGenesWindow = NGenesPopup()
            self.NGenesWindow.show()
            self.NGenesWindow.exec()
            #Thread start
            thread = vizualizationThread(self.adata,self, self.AutoClustering)
            thread.signal.progress_signal.connect(self.reportProgress)
            thread.signal.return_signal.connect(self.completePreprocess)
            self.threadpool.start(thread)
         
    def reportProgress(self,signal):
        print(signal)
        self.parent.parent.parent.statusbar.showMessage(signal)
        self.displayMessageWindow = th.displayMessagePopup(signal)
        self.displayMessageWindow.show()

    def completePreprocess(self):
        self.plotComboBox.setEnabled(True)
        self.displayMessageWindow.close()
        self.parent.parent.parent.statusbar.clearMessage()
        if plot == 'Gene Ranking':
            ax = self.figure.add_subplot(111) 
            img = cv.imread(os.path.join(basedir, 'figures','rank_genes_groups_leiden.png'))
            ax.imshow(img)
        elif  plot == 'Dot Plot':
            ax = self.figure.add_subplot(111) 
            img = cv.imread(os.path.join(basedir, 'figures','dotplot_.png'))
            ax.imshow(img)
        elif plot == 'Violin':
            self.grid.addWidget(self.slider)
            self.updateViolin(0)
        elif plot == 'Stacked Violin':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread(os.path.join(basedir, 'figures','stacked_violin_.png'))
            ax.imshow(img)
        elif plot == 'Matrix Plot':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread(os.path.join(basedir, 'figures','matrixplot_.png'))
            ax.imshow(img)
        elif plot == 'Heatmap':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread(os.path.join(basedir, 'figures','heatmap.png'))
            ax.imshow(img)
        elif plot == 'Tracksplot':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread(os.path.join(basedir, 'figures','tracksplot.png'))
            ax.imshow(img)
        self.canvasPop.draw()

class vizualizationThread(QRunnable):
    signal = pyqtSignal()

    def __init__(self,adata,parent,AutoClustering):
        super(vizualizationThread, self).__init__()     
        self.signal = th.Signals() 
        self.adata = adata
        self.parent = parent
        self.AutoClustering = AutoClustering

    @pyqtSlot()
    def run(self):
        if plot == 'Gene Ranking':
            self.signal.progress_signal.emit("Generating Gene Ranking Graph...")
            sc.pl.rank_genes_groups(self.adata, n_genes=n_genes, sharey=False, key='wilcoxon', show=False, save='.png')
        elif plot == 'Dot Plot':
            self.signal.progress_signal.emit("Generating Dot Plot...")
            sc.pl.rank_genes_groups_dotplot(self.adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', dendrogram = False, show=False, save='.png')
        elif plot == 'Violin':
            self.signal.progress_signal.emit("Generating Violin Graph...")
            sc.pl.rank_genes_groups_violin(self.adata, n_genes=n_genes, key='wilcoxon', show=False, save='.png')
        elif plot == 'Stacked Violin':
            self.signal.progress_signal.emit("Generating Stacked Violin Graph...")
            sc.pl.rank_genes_groups_stacked_violin(self.adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
        elif plot == 'Matrix Plot':
            self.signal.progress_signal.emit("Generating Matrix Plot...")
            sc.pl.rank_genes_groups_matrixplot(self.adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram = False)
        elif plot == 'Heatmap':
            self.signal.progress_signal.emit("Generating Heatmap...")
            if self.AutoClustering:
                sc.pl.rank_genes_groups_heatmap(self.adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show_gene_labels = True, show=False, save='.png', dendrogram = False)
            else:
                adata_sub = self.adata[self.adata.obs['leiden'].isin(['0','1']),:]
                sc.pl.rank_genes_groups_heatmap(adata_sub, n_genes=n_genes, key='wilcoxon', groupby='leiden', show_gene_labels = True, show=False, save='.png', dendrogram = False)
        elif plot == 'Tracksplot':
            self.signal.progress_signal.emit("Generating Tracksplot...")
            if self.AutoClustering:
                sc.pl.rank_genes_groups_tracksplot(self.adata, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram= False)
            else:
                adata_sub = self.adata[self.adata.obs['leiden'].isin(['0','1']),:]
                sc.pl.rank_genes_groups_tracksplot(adata_sub, n_genes=n_genes, key='wilcoxon', groupby='leiden', show=False, save='.png', dendrogram= False)
        self.signal.return_signal.emit()

class NGenesPopup(QDialog):
    def __init__(self):
        QDialog.__init__(self)
        self.setStyleSheet('background-color:white;')
        self.setWindowTitle('NGenes')
        self.textbox = QLineEdit(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
        self.button = QPushButton('OK', self)
        self.button.setStyleSheet('background-color: #7ad0eb;')
        self.button.move(20,80)
        self.button.resize(280,40)
        self.button.clicked.connect(self.on_click)
        grid = QVBoxLayout()
        grid.addWidget(self.textbox)
        grid.addWidget(self.button)

    def on_click(self):
        global n_genes
        if self.textbox.text():
            n_genes = int(self.textbox.text())
        else:
            n_genes = 5
        self.close()