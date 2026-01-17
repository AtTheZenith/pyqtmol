import requests
import json
from PySide6.QtCore import QObject, Signal, QRunnable, Slot

from src.core.molecule import MoleculeService

class MoleculeFetcher(QObject):
    """
    Worker signal interface.
    """
    finished = Signal(str)  # Emits PDB content
    error = Signal(str)     # Emits error message

class MoleculeWorker(QRunnable):
    """
    Background worker to fetch and process molecule data.
    """
    def __init__(self, compound_name: str):
        super().__init__()
        self.compound_name = compound_name
        self.signals = MoleculeFetcher()

    @Slot()
    def run(self):
        try:
            # 1. Check Cache
            cached = MoleculeService.get_cached_structure(self.compound_name)
            if cached:
                self.signals.finished.emit(cached)
                return

            # 2. Try PubChem
            name_safe = self.compound_name.translate(str.maketrans('₀₁₂₃₄₅₆₇₈₉', '0123456789'))
            url_pubchem = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_safe}/SDF"
            res = requests.get(url_pubchem, timeout=10)
            
            if res.status_code == 200:
                # Convert SDF to PDB (CPU heavy)
                pdb = MoleculeService.convert_smiles_to_pdb(res.text)
                MoleculeService.save_to_cache(self.compound_name, pdb)
                self.signals.finished.emit(pdb)
                return

            # 3. Try RCSB PDB
            # First search for ID
            query = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": name_safe},
                },
                "return_type": "entry",
            }
            url_search = f"https://search.rcsb.org/rcsbsearch/v2/query?json={json.dumps(query)}"
            res_search = requests.get(url_search, timeout=10)
            
            if res_search.status_code == 200:
                data = res_search.json()
                results = data.get("result_set", [])
                if results:
                    pdb_id = results[0]["identifier"]
                    url_pdb = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                    res_pdb = requests.get(url_pdb, timeout=10)
                    if res_pdb.status_code == 200:
                        pdb = res_pdb.text
                        MoleculeService.save_to_cache(self.compound_name, pdb)
                        self.signals.finished.emit(pdb)
                        return

            raise ValueError("Compound not found in PubChem or RCSB")

        except Exception as e:
            self.signals.error.emit(str(e))
