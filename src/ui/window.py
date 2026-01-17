import os
from PySide6.QtCore import QThreadPool, Slot
from PySide6.QtGui import QIcon
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QFrame,
    QMessageBox,
    QProgressBar
)

from src.core.config import ASSETS_DIR, CACHE_DIR, ensure_dirs
from src.utils.text import subscript_formula
from src.network.client import MoleculeWorker

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        ensure_dirs()
        self.setWindowTitle("pyqtmol - Production")
        self.resize(1000, 800)
        
        # Setup ThreadPool
        self.thread_pool = QThreadPool()
        self.thread_pool.setMaxThreadCount(2) # Limit concurrency

        # UI Setup
        self._init_ui()
        self._load_styles()

    def _init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)

        # Top Bar
        top_layout = QHBoxLayout()
        self.compound_input = QLineEdit()
        self.compound_input.setPlaceholderText("Enter compound name (e.g. Aspirin)")
        self.compound_input.textChanged.connect(self._handle_text_change)
        self.compound_input.returnPressed.connect(self.search_compound) # Enter key
        
        self.search_button = QPushButton("Search")
        self.search_button.clicked.connect(self.search_compound)
        
        top_layout.addWidget(self.compound_input)
        top_layout.addWidget(self.search_button)
        main_layout.addLayout(top_layout)

        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0) # Indeterminate
        self.progress_bar.hide()
        main_layout.addWidget(self.progress_bar)

        # Web View
        webview_holder = QFrame()
        wv_layout = QVBoxLayout(webview_holder)
        wv_layout.setContentsMargins(0, 0, 0, 0)
        
        self.webview = QWebEngineView()
        # Default empty page
        none_path = os.path.join(ASSETS_DIR, "none.html")
        if os.path.exists(none_path):
             self.webview.setUrl(f"file:///{none_path.replace(os.sep, '/')}")
        
        wv_layout.addWidget(self.webview)
        main_layout.addWidget(webview_holder)

        # Basic Icon
        icon_path = os.path.join(ASSETS_DIR, "64x.ico")
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))

    def _load_styles(self):
        style_path = os.path.join(ASSETS_DIR, "styles.qss")
        if os.path.exists(style_path):
            with open(style_path, "r") as f:
                self.setStyleSheet(f.read())

    def _handle_text_change(self, text):
        formatted = subscript_formula(text)
        if formatted != text:
             self.compound_input.blockSignals(True)
             cursor = self.compound_input.cursorPosition()
             self.compound_input.setText(formatted)
             self.compound_input.setCursorPosition(cursor)
             self.compound_input.blockSignals(False)

    def search_compound(self):
        name = self.compound_input.text().strip()
        if not name:
            return

        # Disable UI
        self.search_button.setEnabled(False)
        self.compound_input.setEnabled(False)
        self.progress_bar.show()

        # Start Worker
        worker = MoleculeWorker(name)
        worker.signals.finished.connect(self._on_search_finished)
        worker.signals.error.connect(self._on_search_error)
        self.thread_pool.start(worker)

    @Slot(str)
    def _on_search_finished(self, pdb_content):
        self._reset_ui()
        self._display_structure(pdb_content)

    @Slot(str)
    def _on_search_error(self, error_msg):
        self._reset_ui()
        QMessageBox.warning(self, "Search Failed", f"Could not find compound:\n{error_msg}")

    def _reset_ui(self):
        self.search_button.setEnabled(True)
        self.compound_input.setEnabled(True)
        self.progress_bar.hide()

    def _display_structure(self, pdb_content):
        # We need to inject this into the template
        # Ideally, we don't write to disk just to read it back, but existing JS logic might rely on it.
        # Original: writes to cache/web.html and loads it.
        # We can do the same to keep 3Dmol.js happy without complex JS injection, 
        # or use setHtml with base url.
        
        # Let's stick to the file approach for stability with local assets (JS bundles).
        template_path = os.path.join(ASSETS_DIR, "template.html")
        if not os.path.exists(template_path):
            QMessageBox.critical(self, "Error", "template.html missing in assets.")
            return

        try:
            with open(template_path, "r") as f:
                template = f.read()
            
            # Escape backslashes for JS string (simple replacement)
            # The original code did: struct_str = "\\n".join(struct.splitlines())
            # This turns newlines into \n literal for JS.
            escaped_pdb = "\\n".join(pdb_content.splitlines())
            html = template.replace("temp", escaped_pdb)
            
            out_path = os.path.join(CACHE_DIR, "web.html")
            with open(out_path, "w") as f:
                f.write(html)
            
            self.webview.setUrl(f"file:///{out_path.replace(os.sep, '/')}")
            
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to render structure: {e}")
