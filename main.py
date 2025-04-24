import json
import os
import re
import sys

import requests
from PySide6.QtGui import QIcon
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import (
    QApplication,
    QFrame,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
from rdkit import Chem
from rdkit.Chem import rdDistGeom

dependencies = {
    "./assets/64x.ico": "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/64x.ico",
    "./assets/none.html": "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/none.html",
    "./assets/styles.qss": "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/styles.qss",
    "./assets/template.html": "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/template.html",
    "./cache/3Dmol-min.js": "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.1/3Dmol-min.js",
}


def check_file(file_path: str, url: str):
    if not os.path.isfile(file_path):
        res = requests.get(url)

        if res.status_code == 200:
            with open(file_path, "wb") as dlf:
                dlf.write(res.content)
                dlf.close()


def subscript_formula(formula: str) -> str:
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    match = re.match(r"^(\d+)?([A-Za-z].*)$", formula)
    if not match:
        return formula

    coeff, compound = match.groups()

    parts = []
    i = 0
    while i < len(compound):
        c = compound[i]
        if c.isalpha() or c in "()":
            parts.append(c)
            i += 1
        elif c.isdigit():
            j = i
            while j < len(compound) and compound[j].isdigit():
                j += 1
            subscript = compound[i:j].translate(subscript_map)
            parts.append(subscript)
            i = j
        else:
            parts.append(c)
            i += 1

    transformed = "".join(parts)
    return (coeff or "") + transformed


def get_file(file: str):
    return str.format("file:///{dir}", dir=os.path.abspath(file).replace("\\", "/"))


def get_struct(compound):
    if os.path.isfile(f"./cache/{compound.upper()}.pdb"):
        with open(f"./cache/{compound.upper()}.pdb", "r") as f:
            return f.read()

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound.translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789'))}/property/CanonicalSMILES/JSON"
    res = requests.get(url)
    if res.status_code == 200:
        try:
            data = res.json()
            smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)

            params = rdDistGeom.ETKDG()
            success = rdDistGeom.EmbedMolecule(mol, params)
            if success != 0:
                raise ValueError("3D embedding failed.")
            pdb = Chem.MolToPDBBlock(mol)
            with open(f"./cache/{compound.upper()}.pdb", "w") as f:
                f.write(pdb)
            return pdb
        except (KeyError, IndexError):
            pass

    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": compound.translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789"))
            },
        },
        "return_type": "entry",
    }

    query_json = json.dumps(query)

    url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={query_json}"
    res = requests.get(url)

    if res.status_code == 200:
        data = res.json()
        pdb_ids = [entry["identifier"] for entry in data.get("result_set", [])]
        pdb_urls = [
            f"https://files.rcsb.org/download/{pdb_id}.pdb" for pdb_id in pdb_ids
        ]

        if len(pdb_ids) != 0:
            pdb = requests.get(pdb_urls[0]).text
            with open(f"./cache/{compound.upper()}.pdb", "w") as f:
                f.write(pdb)
            return pdb

    return None


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("pyqtmol")
        self.showNormal()

        if not os.path.isdir("./assets"):
            os.makedirs("./assets")

        if not os.path.isdir("./cache"):
            os.makedirs("./cache")

        for file_path, url in dependencies.items():
            check_file(file_path, url)

        self.setWindowIcon(QIcon("./assets/64x.ico"))

        main_layout = QVBoxLayout()
        top_layout = QHBoxLayout()
        wv_layout = QVBoxLayout()

        self.compound_input = QLineEdit(self)
        self.compound_input.setPlaceholderText("Enter compound name")
        top_layout.addWidget(self.compound_input)
        self.search_button = QPushButton("Search", self)
        top_layout.addWidget(self.search_button)
        main_layout.addLayout(top_layout)

        webview_holder = QFrame(self)
        self.webview_frame = QWebEngineView(self)
        self.webview_frame.setUrl(get_file("./assets/none.html"))
        wv_layout.addWidget(self.webview_frame)
        webview_holder.setLayout(wv_layout)
        main_layout.addWidget(webview_holder)

        self.compound_input.textChanged.connect(self.subscript_formula)
        self.search_button.clicked.connect(self.search_compound)

        self.setLayout(main_layout)
        with open("./assets/styles.qss", "r") as style:
            self.setStyleSheet(style.read())

    def subscript_formula(self, text):
        formatted = subscript_formula(text)
        if formatted != text:
            cursor_pos = self.compound_input.cursorPosition()
            self.compound_input.blockSignals(True)
            self.compound_input.setText(formatted)
            self.compound_input.setCursorPosition(cursor_pos)
            self.compound_input.blockSignals(False)

    def search_compound(self):
        compound_name = self.compound_input.text().strip()
        struct = get_struct(compound_name)
        if struct:
            self.display_3d_structure(struct)

    def display_3d_structure(self, struct):
        struct_str = "\\n".join(struct.splitlines())

        with open("./assets/template.html", "r") as fetch:
            data = fetch.read()
            data = data.replace("temp", struct_str)
            file = open("./cache/web.html", "w")
            file.write(data)
            file.close()
            fetch.close()

        self.webview_frame.setUrl(get_file("./cache/web.html"))


app = QApplication(sys.argv)
window = MainWindow()
sys.exit(app.exec())
