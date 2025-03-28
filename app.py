import json
import os
from PySide6.QtWidgets import QApplication, QWidget, QVBoxLayout, QLineEdit, QPushButton
from PySide6.QtWebEngineWidgets import QWebEngineView
from rdkit import Chem
from rdkit.Chem import rdDistGeom
import requests
import sys

def get_file(file):
    return str.format("file:///{dir}", dir = os.path.abspath(file).replace("\\", "/"))

def get_struct(compound):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound}/property/CanonicalSMILES/JSON"
    response = requests.get(url)
    
    if not os.path.exists("./cache"):
        os.makedirs("./cache")
    if not os.path.exists("./cache/3Dmol-min.js"):
        with open("./cache/3Dmol-min.js", "w") as f:
            f.write(requests.get("https://3Dmol.org/build/3Dmol-min.js").text)
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

        self.setWindowTitle("Compound SMILES Search")
        self.setGeometry(200, 200, 800, 600)
        
        layout = QVBoxLayout()
        
        self.compound_input = QLineEdit(self)
        self.compound_input.setPlaceholderText("Enter compound name")
        layout.addWidget(self.compound_input)
        
        self.search_button = QPushButton("Search", self)
        layout.addWidget(self.search_button)

        self.web_view = QWebEngineView(self)
        self.web_view.setUrl(get_file('none.html'))
        layout.addWidget(self.web_view)
        
        self.search_button.clicked.connect(self.search_compound)

        self.setLayout(layout)
    
    def search_compound(self):
        compound_name = self.compound_input.text().strip()
        struct = get_struct(compound_name)
        if struct:
            self.display_3d_structure(struct)
    
    def display_3d_structure(self, struct):
        struct_str = '\\n'.join(struct.splitlines())
        
        with open('./template.html', 'r') as fetch:
            template = fetch.read()
            template = template.replace('temp', struct_str)
            file = open('./web.html', 'w')
            file.write(template)
            file.close()
            fetch.close()

        self.web_view.setUrl(get_file('web.html'))


app = QApplication(sys.argv)
window = MainWindow()
window.showFullScreen()
sys.exit(app.exec())
