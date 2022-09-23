import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QLabel, QComboBox, QWidget,QLineEdit
from PyQt5.QtGui import QPainter
from PyQt5.QtCore import QRect
from scipy.io import mmread
import pandas as pd
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.widgets import LassoSelector
from matplotlib.path import Path

class Cluster:
    ClusterName = ""
    ClusterData = []

metadataDF = pd.read_csv("sc_example_data/metadata.csv")
pcaDF = pd.read_csv("sc_example_data/pca.csv")
umapDF = pd.read_csv("sc_example_data/umap.csv")
df = pd.merge(metadataDF,pcaDF, on="Unnamed: 0",  how='inner')
df = pd.merge(df,umapDF, on="Unnamed: 0",  how='inner')
features = pd.read_csv ("sc_example_data/aggr_iHPF_pHPF_N1_SFT/outs/filtered_feature_bc_matrix/features.tsv.gz", sep = '\t', header=None)
barcodes = pd.read_csv ("sc_example_data/aggr_iHPF_pHPF_N1_SFT/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", sep = '\t', header=None)
matrix = mmread("sc_example_data/aggr_iHPF_pHPF_N1_SFT/outs/filtered_feature_bc_matrix/matrix.mtx.gz")
matrixDF = pd.DataFrame.sparse.from_spmatrix(matrix)
matrixDF.columns=barcodes[0]
matrixDF = features.join(matrixDF)
matrixDF = matrixDF.set_index(0)
df = df[df['Unnamed: 0'].isin(matrixDF.head())]
Cluster1 = Cluster()
Cluster2 = Cluster()
selectedArray = []


class SelectFromCollection:

    def __init__(self, ax, collection, x, y, w, gene, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other
        self.x = x
        self.y = y
        self.gene = gene
        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)
        self.w = w

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
            selectedArray.append(df.loc[(df[self.x] == selected[0]) & (df[self.y] == selected[1])]['Unnamed: 0'].values[0])
        for barcode in selectedArray:
            print(self.gene," ",barcode," ",matrixDF.loc[self.gene,barcode])
        self.w = MyPopup()
        self.w.setGeometry(QRect(100, 100, 400, 200))
        self.w.show()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

class MyPopup(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.textbox = QLineEdit(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
        self.textbox.setPlaceholderText("Enter Cluster Name")
        self.button = QPushButton('OK', self)
        self.button.move(20,80)
        # connect button to function on_click
        self.button.clicked.connect(self.on_click)
        grid = QVBoxLayout()
        grid.addWidget(self.textbox)
        grid.addWidget(self.button)

    def on_click(self):
        if not Cluster1.ClusterName:
            Cluster1.ClusterName = self.textbox.text()
            Cluster1.ClusterData = selectedArray
        else :
            Cluster2.ClusterName = self.textbox.text()
            Cluster2.ClusterData = selectedArray
        selectedArray=[]
        self.close()

class Window(QDialog):

    # constructor
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle('Cell Visualization App')
        # self.setFixedSize(QSize(400, 500))
        self.showMaximized()

        
        self.typeComboBox = QComboBox(self)
        self.typeComboBox.addItem("PCA",pcaDF.head())
        self.typeComboBox.addItem("UMAP",umapDF.head())
        self.typeLabel = QLabel("Type:")
        self.typeLabel.setBuddy(self.typeComboBox)

        self.xComboBox = QComboBox(self)
        self.yComboBox = QComboBox()
        self.xComboBox.addItems(self.typeComboBox.itemData(0))
        self.yComboBox.addItems(self.typeComboBox.itemData(0))
        self.typeComboBox.activated.connect(self.clicker)

        self.xLabel = QLabel("X:")
        self.xLabel.setBuddy(self.xComboBox)

        self.yLabel = QLabel("Y:")
        self.yLabel.setBuddy(self.yComboBox)

        self.gComboBox = QComboBox()
        self.gComboBox.addItems(features[0])
        self.gLabel = QLabel("Gene:")
        self.gLabel.setBuddy(self.gComboBox)

        self.c1Label = QLabel(Cluster1.ClusterName)
        self.c2Label = QLabel(Cluster2.ClusterName)
        # self.button = QPushButton('OK', self)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.w = None

        button = QPushButton("Plot Current Attributes", self)
        button.pressed.connect(self.changeValue)

        grid = QVBoxLayout()
        grid.addWidget(self.typeLabel)
        grid.addWidget(self.typeComboBox)
        grid.addWidget(self.xLabel)
        grid.addWidget(self.xComboBox)
        grid.addWidget(self.yLabel)
        grid.addWidget(self.yComboBox)
        grid.addWidget(self.gLabel)
        grid.addWidget(self.gComboBox)
        grid.addWidget(self.c1Label)
        grid.addWidget(self.c2Label)
        grid.addWidget(self.canvas)
        grid.addWidget(button)
        self.setLayout(grid)

    def clicker(self, index):
        self.xComboBox.clear()
        self.yComboBox.clear()
        self.xComboBox.addItems(self.typeComboBox.itemData(index))
        self.yComboBox.addItems(self.typeComboBox.itemData(index))

    def changeValue(self, *args):
        # finding the content of current item in combo box
        plttype = self.typeComboBox.currentText()
        x = self.xComboBox.currentText()
        y = self.yComboBox.currentText()
        gene = self.gComboBox.currentText()

        #clear current figure
        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)

        #create scatter plot with new data
        b = ax.scatter(x=df[x], y=df[y], s = 4)

        #create labels and title
        t = y + " vs " + x
        ax.set(xlabel=x, ylabel =y, title=t )
        selector = SelectFromCollection(ax, b, x, y, self.w, gene)

        #draw new graph
        self.canvas.draw()

if __name__ == '__main__':
        # creating apyqt5 application
        app = QApplication(sys.argv)
        main = Window()
        main.show()
        sys.exit(app.exec_())