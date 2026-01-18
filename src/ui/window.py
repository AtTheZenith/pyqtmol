import os
from PySide6.QtCore import QThreadPool, Slot, Qt
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QMessageBox, QLabel, QCheckBox

from src.core.config import ASSETS_DIR
from src.network.client import SearchWorker, FetchWorker, SearchResult
from src.ui.theme import Theme
from src.ui.widgets import SearchBar, ResultList, MoleculeViewer, LoadingOverlay


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("pyqtmol")
        self.resize(1100, 750)

        self.thread_pool = QThreadPool()
        self.thread_pool.setMaxThreadCount(2)

        Theme.apply(QApplication.instance())
        self._init_ui()

    def _init_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # --- Sidebar ---
        sidebar = QWidget()
        sidebar.setObjectName("sidebar")
        sidebar.setFixedWidth(350)
        # Style moved to QSS
        sidebar_layout = QVBoxLayout(sidebar)
        sidebar_layout.setContentsMargins(15, 20, 15, 20)
        sidebar_layout.setSpacing(15)

        # Title / Brand
        brand_label = QLabel("PyQtMol")
        brand_label.setObjectName("brand_label")
        sidebar_layout.addWidget(brand_label)

        # Search Bar
        self.search_bar = SearchBar()
        self.search_bar.search_triggered.connect(self.search_compound)
        sidebar_layout.addWidget(self.search_bar)

        # Hydrogen Toggle
        self.h_toggle = QCheckBox("Show Hydrogens")
        self.h_toggle.setObjectName("hydrogen_toggle")
        self.h_toggle.setChecked(True)
        self.h_toggle.setCursor(Qt.PointingHandCursor)
        self.h_toggle.stateChanged.connect(self._on_toggle_hydrogens)
        sidebar_layout.addWidget(self.h_toggle)

        # Results Label
        self.results_label = QLabel("Results")
        self.results_label.setObjectName("results_label")
        sidebar_layout.addWidget(self.results_label)

        # Result List
        self.result_list = ResultList()
        self.result_list.itemClicked.connect(self._on_result_clicked)
        sidebar_layout.addWidget(self.result_list)

        # --- Main Content ---
        content_area = QWidget()
        content_layout = QVBoxLayout(content_area)
        content_layout.setContentsMargins(0, 0, 0, 0)

        self.viewer = MoleculeViewer()
        content_layout.addWidget(self.viewer)

        # Splitter (optional if we want resizable sidebar, but fixed is often cleaner for this app)
        # using layout directly for now:
        main_layout.addWidget(sidebar)
        main_layout.addWidget(content_area)

        # --- Overlays ---
        self.loading_overlay = LoadingOverlay(central_widget)

        # Icon
        icon_path = os.path.join(ASSETS_DIR, "icon.ico")
        if os.path.exists(icon_path):
            self.setWindowIcon(QIcon(icon_path))

    def resizeEvent(self, event):
        self.loading_overlay.resize(self.size())
        super().resizeEvent(event)

    def search_compound(self, name):
        self._set_loading(True)
        self.result_list.clear()

        worker = SearchWorker(name)
        worker.signals.finished.connect(self._on_search_finished)
        worker.signals.error.connect(self._on_error)
        self.thread_pool.start(worker)

    def _on_result_clicked(self, item):
        result: SearchResult = item.data(Qt.UserRole)
        if not result:
            return

        self._set_loading(True)

        worker = FetchWorker(result)
        worker.signals.finished.connect(self._on_fetch_finished)
        worker.signals.error.connect(self._on_error)
        self.thread_pool.start(worker)

    @Slot(object)
    def _on_search_finished(self, results):
        self._set_loading(False)
        self.results_label.setText(f"Results ({len(results)})")

        for res in results:
            self.result_list.add_result(res)

    @Slot(str)
    def _on_fetch_finished(self, pdb_content):
        self._set_loading(False)
        self.viewer.display_structure(pdb_content)

    @Slot(str)
    def _on_error(self, error_msg):
        self._set_loading(False)
        QMessageBox.warning(self, "Error", str(error_msg))

    @Slot(int)
    def _on_toggle_hydrogens(self, state):
        self.viewer.set_hydrogens_visible(state == 2)

    def _set_loading(self, loading: bool):
        self.loading_overlay.show_loading(loading)
        self.search_bar.setEnabled(not loading)
        self.result_list.setEnabled(not loading)
