import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QFormLayout, QHBoxLayout, QLabel, QComboBox, QWidget,QLineEdit
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
    def __str__(self): 
        return self.ClusterName+": "+ str(self.ClusterData)

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

    def __init__(self, ax, collection, x, y, w, c1Label, c2Label, alpha_other=0.3):
        self.canvas = ax.figure.canvas
        self.collection = collection
        self.alpha_other = alpha_other
        self.x = x
        self.y = y
        self.xys = collection.get_offsets()
        self.Npts = len(self.xys)
        self.w = w
        self.c1Label = c1Label
        self.c2Label = c2Label

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
        self.w = MyPopup(self.c1Label, self.c2Label)
        self.w.setGeometry(QRect(100, 100, 400, 200))
        self.w.show()

    def disconnect(self):
        self.lasso.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

class MyPopup(QWidget):
    def __init__(self, c1Label, c2Label):
        QWidget.__init__(self)
        self.textbox = QLineEdit(self)
        self.textbox.move(20, 20)
        self.textbox.resize(280,40)
        self.textbox.setPlaceholderText("Enter Cluster Name")
        self.button = QPushButton('OK', self)
        self.button.move(20,80)
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

class Window(QDialog):

    # constructor
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle('Cell Visualization App')
        # self.setFixedSize(QSize(400, 500))
        # self.setStyleSheet("background-color: White;")
        self.showMaximized()

        self.typeComboBox = QComboBox(self)
        self.typeComboBox.addItem("PCA",pcaDF.head())
        self.typeComboBox.addItem("UMAP",umapDF.head())
        self.typeLabel = QLabel("Type:")
        self.typeLabel.setBuddy(self.typeComboBox)
        self.typeComboBox.resize(165, self.typeComboBox.height());

        self.xComboBox = QComboBox(self)
        self.xComboBox.addItems(self.typeComboBox.itemData(0))
        self.xLabel = QLabel("X:")
        self.xLabel.setBuddy(self.xComboBox)

        self.yComboBox = QComboBox()
        self.yComboBox.addItems(self.typeComboBox.itemData(0))
        self.yLabel = QLabel("Y:")
        self.yLabel.setBuddy(self.yComboBox)
        self.typeComboBox.activated.connect(self.clicker)
        
        self.gComboBox = QComboBox()
        self.gComboBox.addItems(features[0])
        self.gLabel = QLabel("Gene:")
        self.gLabel.setBuddy(self.gComboBox)

        plotButton = QPushButton("Plot Current Attributes", self)
        plotButton.pressed.connect(self.changeValue)

        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.w = None

        
        self.c1Label = QLabel(Cluster1.ClusterName)
        self.c2Label = QLabel(Cluster2.ClusterName)
        self.c1Label.setStyleSheet("border: 1px solid black;")
        self.c2Label.setStyleSheet("border: 1px solid black;")
        self.compareButton = QPushButton("Compare", self)
        self.compareButton.pressed.connect(self.compare)
        self.resetButton = QPushButton("Reset", self)
        self.resetButton.pressed.connect(self.reset)

        outerLayout = QVBoxLayout()
        attributeLayout = QFormLayout()
        comparisonLayout = QHBoxLayout()

        attributeLayout.addRow(self.typeLabel,self.typeComboBox)
        attributeLayout.addRow(self.xLabel,self.xComboBox)
        attributeLayout.addRow(self.yLabel,self.yComboBox)
        attributeLayout.addRow(self.gLabel,self.gComboBox)

        comparisonLayout.addWidget(self.c1Label)
        comparisonLayout.addWidget(self.c2Label)
        comparisonLayout.addWidget(self.compareButton)
        comparisonLayout.addWidget(self.resetButton)

        outerLayout.addLayout(attributeLayout)
        outerLayout.addWidget(plotButton)
        outerLayout.addWidget(self.canvas)
        outerLayout.addLayout(comparisonLayout)

        # layout.setContentsMargins(left, top, right, bottom)
        outerLayout.setContentsMargins(50, 20, 50, 50)
        self.setLayout(outerLayout)


    def clicker(self, index):
        self.xComboBox.clear()
        self.yComboBox.clear()
        self.xComboBox.addItems(self.typeComboBox.itemData(index))
        self.yComboBox.addItems(self.typeComboBox.itemData(index))

    def reset(self, *args):
        Cluster1.ClusterName = ""
        Cluster1.ClusterData = []
        Cluster2.ClusterName = ""
        Cluster2.ClusterData = []
        self.c1Label.setText("")
        self.c2Label.setText("")
        self.changeValue()

    def compare(self, *args):
        gene = self.gComboBox.currentText()
        self.figure.clear() #clear current figure
        ax = self.figure.add_subplot(121) # create an axis
        c1matrix = []
        c2matrix = []
        for barcode in Cluster1.ClusterData:
            c1matrix.append(matrixDF.loc[gene,barcode])
            print(Cluster1.ClusterName,": ",gene," ",barcode," ",matrixDF.loc[gene,barcode])
        ax.hist(c1matrix)
        ax = self.figure.add_subplot(122) # create an axis
        for barcode in Cluster2.ClusterData:
            c2matrix.append(matrixDF.loc[gene,barcode])
            print(Cluster2.ClusterName,": ",gene," ",barcode," ",matrixDF.loc[gene,barcode])
        ax.hist(c2matrix)
        self.canvas.draw()


    def changeValue(self, *args):

        Cluster1.ClusterName = ""
        Cluster1.ClusterData = []
        Cluster2.ClusterName = ""
        Cluster2.ClusterData = []
        self.c1Label.setText("")
        self.c2Label.setText("")
        
        # finding the content of current item in combo box
        plttype = self.typeComboBox.currentText()
        x = self.xComboBox.currentText()
        y = self.yComboBox.currentText()
        self.figure.clear() #clear current figure
        ax = self.figure.add_subplot(111) # create an axis
        b = ax.scatter(x=df[x], y=df[y], s = 4) #create scatter plot with new data

        #create labels and title
        t = y + " vs " + x
        ax.set(xlabel=x, ylabel =y, title=t )
        selector = SelectFromCollection(ax, b, x, y, self.w, self.c1Label, self.c2Label)

        #draw new graph
        self.canvas.draw()

if __name__ == '__main__':
        # creating apyqt5 application
        app = QApplication(sys.argv)
        main = Window()
        main.show()
        sys.exit(app.exec_())