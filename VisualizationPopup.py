import cv2 as cv
import matplotlib.pyplot as plt
import scanpy as sc
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import QComboBox, QLabel, QVBoxLayout, QWidget, QFormLayout,QPushButton,QHBoxLayout

import math


class visualizationPopup(QWidget):
    def __init__(self,grouplen,AutoClustering):
        QWidget.__init__(self)
        self.grouplen=grouplen
        self.setWindowTitle('Differential Gene Expression Representation')
        self.showMaximized()
        self.plotComboBox = QComboBox()
        if AutoClustering:
            self.plotComboBox.addItems(['Gene Ranking','Dendrogram', 'Dot Plot', 'Violin', 'Stacked Violin', 'Matrix Plot', 'Heatmap', 'Tracksplot'])
        else:
            self.plotComboBox.addItems(['Gene Ranking','Dendrogram', 'Dot Plot', 'Violin', 'Stacked Violin', 'Matrix Plot'])
        self.plotLabel = QLabel("Plot Type:")
        self.plotLabel.setBuddy(self.plotComboBox)
        
        self.plotComboBox.currentIndexChanged.connect(self.updateGraph)

        #ELEPHANT Make rank grapg default
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)

        self.toolbar = NavigationToolbar(self.canvas, self)
        self.plotSelectionLayout = QFormLayout()
        self.plotSelectionLayout.addRow(self.plotLabel,self.plotComboBox)

        grid = QVBoxLayout()
        grid.addLayout(self.plotSelectionLayout)
        grid.addWidget(self.toolbar)
        grid.addWidget(self.canvas,1)
        self.setLayout(grid)
        

    def updateGraph(self, *args):

        plot = self.plotComboBox.currentText()
        self.figure.clear() #clear current figure
        self.canvas.flush_events()
        if plot == 'Gene Ranking':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/rank_genes_groups_leiden.png')
            ax.imshow(img)
        elif plot == 'Dendrogram':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/dendrogram.png')
            ax.imshow(img)
        elif  plot == 'Dot Plot':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/dotplot_.png')
            ax.imshow(img)
        #ELEPHANT figure out the plots
        elif plot == 'Violin':
            i=0
            index=1
            if(self.grouplen<=3):
                column = self.grouplen
                row=1
            else:
                column = 3
                row = math.ceil(self.grouplen/3)
            
            for i in range (self.grouplen):
                ax = self.figure.add_subplot(row,column,index)
                index=index+1
                img = cv.imread('./figures/rank_genes_groups_leiden_'+str(i)+'.png')
                ax.imshow(img)
        elif plot == 'Stacked Violin':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/stacked_violin_.png')
            ax.imshow(img)
        elif plot == 'Matrix Plot':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/matrixplot_.png')
            ax.imshow(img)
        elif plot == 'Heatmap':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/heatmap.png')
            ax.imshow(img)
        elif plot == 'Tracksplot':
            ax = self.figure.add_subplot(111) # create an axis
            img = cv.imread('./figures/tracksplot.png')
            ax.imshow(img)
        self.canvas.draw()