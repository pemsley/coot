import libBoostedMoleculesToTrianglesCXXClasses as CXXClasses
from libBoostedMoleculesToTrianglesCXXClasses import SceneSetup, Camera, CXXCoord_float, CompoundSelection, ColorScheme, MolecularRepresentationInstance, MyMolecule, Light, RendererGLSL, RendererGL, SolidColorRule
import sys

import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QPushButton, QMainWindow, QTextEdit
from PyQt5.QtOpenGL import QGLWidget
from PyQt5.QtGui import (QGuiApplication, QMatrix4x4, QOpenGLContext,
                         QOpenGLShader, QOpenGLShaderProgram, QSurfaceFormat, QWindow, QMouseEvent)
from PyQt5.QtCore import QPoint, Qt
from PyQt5.QtCore import QObject, pyqtSignal
from PyQt5.QtWidgets import QWidget, QGridLayout, QLineEdit, QLabel
import math

class WfWidget(QGLWidget):
    def __init__(self, parent = None):
        super(WfWidget, self).__init__(parent)
        self.prepareCXXClasses()
        self.setMouseTracking(True)
        self.oldMousePosition = QPoint(-1,-1)
        self.tested = False
    
    def mouseMoveEvent(self, event):
        newMousePosition = event.pos()
        deltaX = float(newMousePosition.x() - self.oldMousePosition.x())
        deltaY = float(newMousePosition.y() - self.oldMousePosition.y())
        magnitude = math.sqrt((deltaX*deltaX) + (deltaY*deltaY))
        if magnitude > 0:
            if event.buttons() & Qt.LeftButton :
                if event.modifiers() == Qt.NoModifier:
                    xComponent =  deltaY / magnitude
                    yComponent =  deltaX / magnitude
                    rotation = CXXCoord_float(magnitude/4., xComponent, yComponent, 0.)
                    self.camera.getSceneSetup().rotateBy(rotation)
                    self.updateGL()
            elif event.buttons() & Qt.RightButton :
                oldFovy = self.camera.getFovy()
                factor = 1. + (abs(deltaY)/100.)
                newFovy = oldFovy * factor
                if deltaY < 0: newFovy = oldFovy / factor
                self.camera.setFovy(newFovy)
                self.updateGL()
    
        self.oldMousePosition = newMousePosition
    
    def paintGL(self):
        self.makeCurrent()
        self.camera.renderWithRenderer(self.renderer)
    
    def resizeGL(self, w, h):
        pass

    def initializeGL(self):
        super(WfWidget,self).initializeGL()
        self.makeCurrent()
        print(self.context())
        self.makeCurrent()
        self.renderer = RendererGLSL()
        self.renderer.init()
        
        #self.renderer = RendererGL()
    
    def prepareCXXClasses(self):
        light = CXXClasses.Light.defaultLight()
        light.setDrawLight(False)
        sceneSetup = CXXClasses.SceneSetup.defaultSceneSetup()
        sceneSetup.addLight(light)
        self.camera = CXXClasses.Camera.defaultCamera()
        sceneSetup.addCamera(self.camera)
        self.camera.setSceneSetup(sceneSetup)

import sys
import code
from io import StringIO

class Interpreter(QObject, code.InteractiveConsole):
    output = pyqtSignal(str)

    def __init__(self, locals=None):
        QObject.__init__(self)
        code.InteractiveConsole.__init__(self, locals)
        self.out = StringIO()

    def write(self, data):
        self.output.emit(data)

    def runcode(self, codez):
        """
        Reimplementation to capture stdout and stderr
        """
        sys.stdout = self.out
        sys.stderr = self.out
        result = code.InteractiveConsole.runcode(self, codez)
        
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        sys.stdout.flush()
        sys.stderr.flush()
        self.output.emit(self.out.getvalue())

class MainWidget(QWidget):
    from PyQt5 import QtCore
    def __init__(self, *args, **kws):
        super(MainWidget,self).__init__(*args, **kws)
        layout = QGridLayout(self)
        self.display = QTextEdit()
        self.display.setReadOnly(True)
        layout.addWidget(self.display, 2,0, 1,2)
        self.display.setMaximumHeight(100);
        self.history = QTextEdit()
        self.history.setReadOnly(True)
        layout.addWidget(self.history, 2,2, 1,2)
        self.history.setMaximumHeight(100);
        self.prompt = QLabel(">>>")
        self.prompt.setMaximumWidth(100);
        layout.addWidget(self.prompt, 1,0)
        self.input = QLineEdit()
        self.input.installEventFilter(self)
        layout.addWidget(self.input, 1,1, 1, 3)
        self.wfWidget = WfWidget()
        self.wfWidget.setMinimumHeight(300);
        layout.addWidget(self.wfWidget, 0,0, 1,4)
        rotX=CXXCoord_float(1.,1.,0.,0.)
        myDict = {"camera":self.wfWidget.camera,"sceneSetup":self.wfWidget.camera.getSceneSetup()}
        combinedDict = {**globals(), **locals(), **myDict}
        self.interp = Interpreter(locals=combinedDict)
        self.input.returnPressed.connect(self.text_input)
        self.interp.output.connect(self.text_output)
        self.interp.runcode('''exec(open("QMGPy.py").read())''')
        self.inputLines = []
        self.inputLinesIndex = 0

    def eventFilter(self, object, event):
        from PyQt5 import QtCore
        if object == self.input:
            if event.type() == QtCore.QEvent.KeyPress:
                if event.key() == QtCore.Qt.Key_Up:
                    #print( "lineEdit -> Qt::Key_Up")
                    self.inputLinesIndex -= 1
                    if -1*self.inputLinesIndex <= len(self.inputLines):
                        replacementText = self.inputLines[self.inputLinesIndex]
                        self.input.setText(replacementText)
                    else:
                        self.inputLinesIndex = -1*len(self.inputLines)
                    return True
                elif event.key() == QtCore.Qt.Key_Down:
                    #print( "lineEdit -> Qt::Key_Down")
                    self.inputLinesIndex += 1
                    if self.inputLinesIndex <= -1:
                        replacementText = self.inputLines[self.inputLinesIndex]
                        self.input.setText(replacementText)
                    else:
                        self.inputLinesIndex = -1
                    return True
            return False
        return super(self, QWidget).eventFilter(obj, event)

    def text_input(self):
        text = self.input.text()
        self.inputLines.append(text)
        self.inputLinesIndex = 0
        self.history.append(">" + str(text))
        self.history.verticalScrollBar().setValue(self.history.verticalScrollBar().maximum());
        self.input.clear()
        
        if self.interp.push(str(text)):
            # More input required
            # Use sys.ps1 and sys.ps2
            self.prompt.setText("...")
        else:
            self.prompt.setText(">>>")
        self.wfWidget.update()


    def text_output(self, text):
        print(text)
        sys.stdout.flush()
        self.display.setPlainText(text)
        #self.display.append(text)
        self.display.verticalScrollBar().setValue(self.display.verticalScrollBar().maximum());


from PyQt5.QtWidgets import QTextEdit, QPushButton, QVBoxLayout

if __name__ == '__main__':
    app = QApplication(["CXXClasses"])
    mainWindow = MainWidget()
    mainWindow.show()
    sys.exit(app.exec_())
