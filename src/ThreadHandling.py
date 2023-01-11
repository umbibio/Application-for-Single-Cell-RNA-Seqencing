from PyQt5.QtGui import QMovie
from PyQt5.QtCore import QObject,pyqtSignal,QRect,Qt
from PyQt5.QtWidgets import QDialog, QLabel,QHBoxLayout

class Signals(QObject):
    return_signal = pyqtSignal()
    progress_signal = pyqtSignal(str)

class displayMessagePopup(QDialog):
    def __init__(self,signal):
        QDialog.__init__(self)
        self.setWindowFlag(Qt.FramelessWindowHint)
        # self.setStyleSheet('background-color:light white;')
        # self.setWindowOpacity(0.7)
        self.setStyleSheet('background-color:light grey;')
        self.message = QLabel(self)
        self.message.setText(signal)
        self.message.move(60, 25)
        self.message.resize(500,20)
        self.message.setStyleSheet("border: None")
        self.MovieLabel = QLabel(self)
        self.MovieLabel.setGeometry(QRect(0, 0, 30, 30))
        self.movie = QMovie("images/loader.gif")
        self.MovieLabel.setMovie(self.movie)
        self.movie.start()
        self.MovieLabel.move(20, 20)
        
        grid = QHBoxLayout()
        grid.addWidget(self.MovieLabel,1) 
        grid.addWidget(self.message,1) 
        self.startAnimation()
  
    # Start Animation
    def startAnimation(self):
        self.movie.start()