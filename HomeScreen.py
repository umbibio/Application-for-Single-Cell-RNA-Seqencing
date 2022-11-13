import MainWindowFunctions as mwf
import DifferentialGeneAnalysis as dgf
import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout,QFileDialog
import scanpy as sc
import os

class Window(QDialog):

    # constructor
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle('Cell Visualization App')
        # self.setFixedSize(QSize(400, 500))
        # self.setStyleSheet("background-color: White;")
        self.showMaximized()

        self.uploadButton = QPushButton("Select folder consisting of .mtx and .tsv files", self)
        self.uploadButton.pressed.connect(self.upload)
        self.uploadh5adButton = QPushButton("Upload new h5ad files", self)
        self.uploadh5adButton.pressed.connect(self.uploadh5ad)
        self.existingFileButton = QPushButton("Use existing file", self)
        self.existingFileButton.pressed.connect(self.existingFile)

        self.outerLayout = QVBoxLayout()
        self.outerLayout.addWidget(self.uploadButton)
        self.outerLayout.addWidget(self.uploadh5adButton)
        self.outerLayout.addWidget(self.existingFileButton)
        self.outerLayout.setContentsMargins(50, 20, 50, 50)
        self.setLayout(self.outerLayout)

    def upload(self, *args):
        try:
            folderpath  = QFileDialog.getExistingDirectory(self,'Select folder consisting of .mtx and .tsv files')
            adata = sc.read_10x_mtx(
                    folderpath,  # the directory with the `.mtx` file
                    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
                    cache=False) 
            adata = mwf.preprocessAnnData(adata)
            self.mainWindow = dgf.Window(adata)
            self.mainWindow.show()
        except Exception as e:
             print('ERROR: ',e)
    
    
    def uploadh5ad(self, *args):
        try:
            filename = QFileDialog.getOpenFileName()
            path = filename[0]
            adata = sc.read(path)
            adata = mwf.preprocessAnnData(adata)
            self.mainWindow = dgf.Window(adata)
            self.mainWindow.show()
        except Exception as e:
             print('ERROR: ',e)


    def existingFile(self, *args):
        try:
            if os.path.isfile('./PreprocessedData/adata.h5ad'):
                adata = sc.read('./PreprocessedData/adata.h5ad') 
                adata.uns['log1p']["base"] = None
            else:
                adata = sc.read_10x_mtx(
                    'sc_example_data/aggr_iHPF_pHPF_N1_SFT/outs/filtered_feature_bc_matrix/',  # the directory with the `.mtx` file
                    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
                    cache=False) 
                adata = mwf.preprocessAnnData(adata)
            self.mainWindow = dgf.Window(adata)
            self.mainWindow.show()
        except Exception as e:
            print('ERROR: ',e)

if __name__ == '__main__':
    # creating apyqt5 application
    app = QApplication(sys.argv)
    # app.setStyleSheet("QLabel{font-size: 1pt;}")
    main = Window()
    main.show()
    sys.exit(app.exec_())