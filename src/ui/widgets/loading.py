from PySide6.QtWidgets import QWidget, QLabel, QVBoxLayout
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QPalette


class LoadingOverlay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_TransparentForMouseEvents)  # Allow clicks through if hidden, but we usually show/hide
        self.hide()

        # Dim background
        self.setAutoFillBackground(True)
        palette = self.palette()
        palette.setColor(QPalette.Window, QColor(0, 0, 0, 150))
        self.setPalette(palette)

        layout = QVBoxLayout(self)
        layout.setAlignment(Qt.AlignCenter)

        self.label = QLabel("Loading...")
        self.label.setObjectName("loading_label")
        layout.addWidget(self.label)

    def show_loading(self, show=True):
        if show:
            self.show()
            self.raise_()
        else:
            self.hide()

    def resizeEvent(self, event):
        self.resize(self.parent().size())
        super().resizeEvent(event)
