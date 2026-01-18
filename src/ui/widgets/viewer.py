import os
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QFrame, QVBoxLayout, QMessageBox
from PySide6.QtCore import QUrl, Qt

from src.core.config import TEMPLATES_DIR, CACHE_DIR


class MoleculeViewer(QFrame):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)

        self.webview = QWebEngineView()
        self.webview.page().setBackgroundColor(Qt.transparent)

        self.layout.addWidget(self.webview)

        self._show_hydrogens = True  # Default state
        self._current_pdb_content = None
        self.load_empty()

    def set_hydrogens_visible(self, visible: bool):
        self._show_hydrogens = visible
        js = f"toggleHydrogens({str(visible).lower()});"
        self.webview.page().runJavaScript(js)

    def load_empty(self):
        none_path = os.path.join(TEMPLATES_DIR, "none.html")
        if os.path.exists(none_path):
            self.webview.setUrl(QUrl.fromLocalFile(none_path))

    def display_structure(self, pdb_content):
        self._current_pdb_content = pdb_content
        template_path = os.path.join(TEMPLATES_DIR, "template.html")
        if not os.path.exists(template_path):
            QMessageBox.critical(self, "Error", "template.html missing")
            return

        try:
            with open(template_path, "r") as f:
                template = f.read()

            # Basic escaping for JS string
            # Use JSON dumps to safely escape the string for JS
            import json

            pdb_json = json.dumps(pdb_content)

            # We inject the JSON string directly into the JS variable declaration
            # The template has: let currentPDB = PDB_CONTENT_PLACEHOLDER;
            html = template.replace("PDB_CONTENT_PLACEHOLDER", pdb_json)

            # Inject Hydrogen State
            # template has: let showHydrogens = HYDROGEN_STATE_PLACEHOLDER;
            html = html.replace("HYDROGEN_STATE_PLACEHOLDER", "true" if self._show_hydrogens else "false")

            out_path = os.path.join(CACHE_DIR, "web.html")
            with open(out_path, "w") as f:
                f.write(html)

            # If we are already viewing this file, reload it to pick up changes on disk
            current_url = self.webview.url()
            new_url = QUrl.fromLocalFile(out_path)

            if current_url.toLocalFile() == new_url.toLocalFile():
                self.webview.reload()
            else:
                self.webview.setUrl(new_url)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to render: {e}")
