from PySide6.QtWidgets import QLineEdit, QFrame, QHBoxLayout, QPushButton
from PySide6.QtCore import Qt, Signal


class SearchBar(QFrame):
    search_triggered = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFixedHeight(50)
        self.setObjectName("search_bar")

        layout = QHBoxLayout(self)
        layout.setContentsMargins(20, 5, 5, 5)

        self.input = QLineEdit()
        self.input.setObjectName("search_input")
        self.input.setPlaceholderText("Search PubChem (e.g. Aspirin)...")
        self.input.returnPressed.connect(self._on_search)

        self.search_btn = QPushButton("Go")
        self.search_btn.setObjectName("search_button")
        self.search_btn.setCursor(Qt.PointingHandCursor)
        self.search_btn.setFixedSize(40, 40)
        self.search_btn.clicked.connect(self._on_search)

        layout.addWidget(self.input)
        layout.addWidget(self.search_btn)

    def _on_search(self):
        text = self.input.text().strip()
        if text:
            self.search_triggered.emit(text)

    def text(self):
        return self.input.text()

    def setText(self, text):
        self.input.setText(text)

    def setEnabled(self, enabled):
        self.input.setEnabled(enabled)
        self.search_btn.setEnabled(enabled)
