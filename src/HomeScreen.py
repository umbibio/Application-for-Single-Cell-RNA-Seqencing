import os
import sys
import scanpy as sc
import anndata as anndata
import ThreadHandling as th
import DifferentialGeneAnalysis as dgf
from PyQt5.QtCore import QFile, QTextStream,Qt, pyqtSignal, QThreadPool, pyqtSlot, QRunnable
from PyQt5.QtWidgets import QApplication,QFileDialog,QMenu, QAction,QMainWindow,QWidget,QLabel,QFrame,QVBoxLayout,QPushButton,QComboBox,QLineEdit,QDesktopWidget

basedir = os.path.dirname(__file__)

adata = anndata.AnnData(
  X = None,
  obs = None,
  var = None)
svd_solver=''
n_neighbors=15

class PreProcessingThread(QRunnable):
    signal = pyqtSignal()

    def __init__(self,parent):
        super(PreProcessingThread, self).__init__()     
        self.signal = th.Signals()  
        self.parent = parent  

    @pyqtSlot()
    def run(self):
        global adata
        self.signal.progress_signal.emit("Preprocessing annData file...")
        # adata = mwf.preprocessAnnData(adata,svd_solver,n_neighbors)
        adata.var_names_make_unique()
        self.signal.progress_signal.emit("Normalize annData file...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        self.signal.progress_signal.emit("Logarithmize annData file...")
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        self.signal.progress_signal.emit("Perform principle component analysis...")
        sc.tl.pca(adata, svd_solver=svd_solver)
        self.signal.progress_signal.emit("Compute neighbors...")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, knn=True, n_pcs=40)
        self.signal.progress_signal.emit("Uniform Manifold Approximation and Projection for Dimension Reduction...")
        sc.tl.umap(adata)
        self.signal.return_signal.emit()

class PreProcessingPopup(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        self.parent=MainWindow()
        self.displayMessageWindow = None
        self.center()
        self.setStyleSheet('background-color:white;')
        self.setWindowTitle('PreProcessing Specifications')
        self.pcaLabel = QLabel("Principal component analysis",self)
        self.pcaLabel.move(20, 20)
        self.pcaLabel.resize(280,40)
        self.SVDLabel = QLabel("SVD solver to use:", self)
        self.SVDLabel.move(20, 80)
        self.SVDLabel.resize(280,20)
        self.svd_solver = QComboBox(self)
        self.svd_solver.addItems(['arpack','randomized', 'auto', 'lobpcg'])
        self.SVDLabel.setBuddy(self.svd_solver)
        self.svd_solver.move(20, 100)
        self.svd_solver.resize(280,40)
        
        self.Separator = QFrame(self)
        self.Separator.setFrameShape(QFrame.HLine)
        self.Separator.setLineWidth(1)
        self.Separator.move(20, 180)
        self.Separator.resize(500,40)

        self.umapLabel = QLabel("Uniform Manifold Approximation and Projection for Dimension Reduction",self)
        self.umapLabel.move(20, 240)
        self.umapLabel.resize(500,40)
        self.neighborsLabel = QLabel("Number of Neighbors:",self)
        self.neighborsLabel.move(20, 300)
        self.neighborsLabel.resize(280,20)
        self.n_neighbors = QLineEdit(self)
        self.n_neighbors.move(20, 320)
        self.n_neighbors.resize(280,40)
        
        self.okButton = QPushButton('OK', self)
        self.okButton.setStyleSheet('background-color: #7ad0eb;')
        self.okButton.move(20,420)
        self.okButton.resize(280,40)
        self.okButton.clicked.connect(self.on_click)
        grid = QVBoxLayout()
        grid.addWidget(self.pcaLabel)
        grid.addWidget(self.SVDLabel)
        grid.addWidget(self.svd_solver)
        grid.addWidget(self.Separator)
        grid.addWidget(self.umapLabel)
        grid.addWidget(self.neighborsLabel)
        grid.addWidget(self.n_neighbors)
        grid.addWidget(self.okButton)
        self.threadpool = QThreadPool()
    
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    def completePreprocess(self):
        self.centralwidget = dgf.Window(self,adata)
        self.parent.setCentralWidget(self.centralwidget)
        self.parent.showMaximized()
        self.parent.show()
        self.displayMessageWindow.close()
        self.close()
    
    def reportProgress(self,signal):
        print(signal)
        self.displayMessageWindow = th.displayMessagePopup(signal)
        self.displayMessageWindow.show()

    def on_click(self):
        self.okButton.setEnabled(False)
        global svd_solver,n_neighbors,adata
        svd_solver = self.svd_solver.currentText()
        n_neighbors = int(self.n_neighbors.text())
        #Thread start
        thread = PreProcessingThread(self)
        thread.signal.progress_signal.connect(self.reportProgress)
        thread.signal.return_signal.connect(self.completePreprocess)
        self.threadpool.start(thread)
             

class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle('Cell Visualization App')
        self.setWindowFlag(Qt.WindowMinimizeButtonHint, True)
        self.setWindowFlag(Qt.WindowMaximizeButtonHint, True)
        self._createActions()
        self._createMenuBar()
        self._connectActions()
        self._createStatusBar()
        self.PreProcessingWindow = None

    def upload(self, *args):
        global adata
        folderpath  = QFileDialog.getExistingDirectory(self,'Select folder consisting of mtx and .tsv files')
        try:
            adata = sc.read_10x_mtx(
                    folderpath,
                    var_names='gene_symbols',
                    cache=False) 
            self.hide()
            self.PreProcessingWindow = PreProcessingPopup()
            # self.PreProcessingWindow.setGeometry(610, 250, 700, 500)
            self.PreProcessingWindow.show()
        except Exception as e:
             print('ERROR: ',e)

    def uploadh5ad(self):
        global adata
        try:
            filename = QFileDialog.getOpenFileName()
            path = filename[0]
            adata = sc.read(path)
            self.hide()
            self.PreProcessingWindow = PreProcessingPopup()
            # self.PreProcessingWindow.setGeometry(610, 250, 700, 500)
            self.PreProcessingWindow.show()
        except Exception as e:
             print('ERROR: ',e)
    
    def _createMenuBar(self):
        menuBar = self.menuBar()
        fileMenu = QMenu("&File", self)
        menuBar.addMenu(fileMenu)
        fileMenu.addAction(self.openh5adAction)
        fileMenu.addAction(self.openFolderAction)
        # fileMenu.addAction(self.saveAction)
        fileMenu.addAction(self.exitAction)
        editMenu = menuBar.addMenu("&Edit")
        helpMenu = menuBar.addMenu("&Help")
    
    def _createStatusBar(self):
        self.statusbar = self.statusBar()

    def _createActions(self):
        self.openh5adAction = QAction(self)
        self.openh5adAction.setText("&Open new h5ad file")
        self.openFolderAction = QAction("&Open folder with mtx and tsv files", self)
        # self.saveAction = QAction("&Save", self)
        self.exitAction = QAction("&Exit", self)
    
    def _connectActions(self):
        self.openh5adAction.triggered.connect(self.uploadh5ad)
        self.openFolderAction.triggered.connect(self.upload)
        # self.saveAction.triggered.connect(self.saveFile)
        self.exitAction.triggered.connect(self.close)
    
    def close(self):
        sys.exit(app.exec_())

if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyleSheet("QLabel{font-size: 1pt;}")
    f = QFile('style.qss')                                
    f.open(QFile.ReadOnly | QFile.Text)
    ts = QTextStream(f)
    stylesheet = ts.readAll()    
    app.setStyleSheet(stylesheet)
    if os.path.isfile(os.path.join(basedir,'PreprocessedData','adata.h5ad')):
        adata = sc.read(os.path.join(basedir,'PreprocessedData','adata.h5ad')) 
        adata.uns['log1p']["base"] = None
    else:
        filename = QFileDialog.getOpenFileName()
        path = filename[0]
        adata = sc.read(path)
    PreProcessingWindow = PreProcessingPopup()
    # PreProcessingWindow.setGeometry(610, 250, 700, 500)
    PreProcessingWindow.show()
    sys.exit(app.exec_())


