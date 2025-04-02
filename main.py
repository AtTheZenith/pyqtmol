import json
import os
from PySide6.QtGui import QIcon
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWidgets import QApplication, QFrame, QHBoxLayout, QLineEdit, QPushButton, QVBoxLayout, QWidget
from rdkit import Chem
from rdkit.Chem import rdDistGeom
import requests
import sys

def get_file(file: str):
    return str.format("file:///{dir}", dir = os.path.abspath(file).replace("\\", "/"))

def get_struct(compound):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if not os.path.exists("./cache"):
        os.makedirs("./cache")
    if not os.path.exists("./cache/3Dmol-min.js"):
        with open("./cache/3Dmol-min.js", "w") as f:
            f.write(requests.get("https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.1/3Dmol-min.js").text)
    if os.path.exists(f"./cache/{compound.lower()}.pdb"):
        with open(f"./cache/{compound.lower()}.pdb", "r") as f:
            return f.read()
    
    if response.status_code == 200:
        try:
            data = response.json()
            smiles = data["PropertyTable"]["Properties"][0]["CanonicalSMILES"]
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol)
            
            params = rdDistGeom.ETKDG()
            success = rdDistGeom.EmbedMolecule(mol, params)
            if success != 0:
                raise ValueError('3D embedding failed.')
            pdb = Chem.MolToPDBBlock(mol)
            with open(f"./cache/{compound.lower()}.pdb", "w") as f:
                f.write(pdb)
            return pdb
        except (KeyError, IndexError):
            pass

    query = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": compound
            }
        },
        "return_type": "entry"
    }

    query_json = json.dumps(query)

    url = f"https://search.rcsb.org/rcsbsearch/v2/query?json={query_json}"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        pdb_ids = [entry['identifier'] for entry in data.get('result_set', [])]
        pdb_urls = [f"https://files.rcsb.org/download/{pdb_id}.pdb" for pdb_id in pdb_ids]
        
        if len(pdb_ids) != 0:
            pdb = requests.get(pdb_urls[1]).text
            with open(f"./cache/{compound.lower()}.pdb", "w") as f:
                f.write(pdb)
            return pdb

    return None


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("pyqtmol")
        self.setWindowIcon(QIcon("./assets/64x.ico"))
        self.showNormal()
        
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
        self.webview_frame.setUrl(get_file('./assets/none.html'))
        wv_layout.addWidget(self.webview_frame)
        webview_holder.setLayout(wv_layout)
        main_layout.addWidget(webview_holder)
        
        self.search_button.clicked.connect(self.search_compound)

        self.setLayout(main_layout)
        with open("./assets/styles.qss", "r") as style:
            self.setStyleSheet(style.read())
    
    def search_compound(self):
        compound_name = self.compound_input.text().strip()
        struct = get_struct(compound_name)
        if struct:
            self.display_3d_structure(struct)
    
    def display_3d_structure(self, struct):
        struct_str = '\\n'.join(struct.splitlines())
        
        with open('./assets/template.html', 'r') as fetch:
            data = fetch.read()
            data = data.replace('temp', struct_str)
            file = open('./cache/web.html', 'w')
            file.write(data)
            file.close()
            fetch.close()

        self.webview_frame.setUrl(get_file('./cache/web.html'))


app = QApplication(sys.argv)
window = MainWindow()
sys.exit(app.exec())