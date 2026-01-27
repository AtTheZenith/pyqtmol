import requests
import json
from dataclasses import dataclass
from typing import List
from PySide6.QtCore import QObject, Signal, QRunnable, Slot

from src.core.molecule import MoleculeService


@dataclass
class SearchResult:
    identifier: str
    name: str
    source: str  # 'PubChem' or 'PDB' or 'Cache'
    inquiry_name: str = ""  # The search term used to find this result


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.
    """

    finished = Signal(object)  # Can be list[SearchResult] or str (pdb content)
    error = Signal(str)


class SearchWorker(QRunnable):
    """
    Searches for compounds and returns a list of matches.
    """

    def __init__(self, query: str):
        super().__init__()
        self.query = query
        self.signals = WorkerSignals()

    @Slot()
    def run(self):
        results: List[SearchResult] = []
        try:
            name_safe = self.query.strip()

            # 1. Check Cache (Simple check for exact match file)
            cached_pdb = MoleculeService.get_cached_structure(name_safe)
            if cached_pdb:
                results.append(SearchResult(identifier=name_safe, name=f"{name_safe} (Cached)", source="Cache", inquiry_name=name_safe))

            # 2. PubChem "Search" (Check if the name itself is valid)
            url_pc = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_safe}/description/JSON"
            res_pc = requests.get(url_pc, timeout=5)
            if res_pc.status_code == 200:
                # If we get a valid response, it exists.
                # PubChem usually resolves 'Aspirin' to a specific CID.
                # We'll just add the query name as a valid result.
                results.append(SearchResult(identifier=name_safe, name=name_safe.capitalize(), source="PubChem", inquiry_name=name_safe))

            # 3. Search RCSB PDB
            # We search PDB first for a list because PubChem often just returns one "best" CID via the name endpoint
            # To get a *list* from PubChem is harder with just a name.
            # So we will do:
            # - PDB Search (lists entries)
            # - PubChem (if name matches exactly-ish)

            # RCSB Search
            q_json = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {"value": name_safe},
                },
                "return_type": "entry",
            }
            url_search = f"https://search.rcsb.org/rcsbsearch/v2/query?json={json.dumps(q_json)}"
            res_search = requests.get(url_search, timeout=5)
            if res_search.status_code == 200:
                data = res_search.json()
                for entry in data.get("result_set", [])[:10]:  # Limit to 10
                    results.append(SearchResult(identifier=entry["identifier"], name=entry["identifier"], source="PDB", inquiry_name=name_safe))

        except Exception as e:
            self.signals.error.emit(str(e))

        if not results or len(results) == 0:
            raise ValueError("No results found.")

        self.signals.finished.emit(results)
            


class FetchWorker(QRunnable):
    """
    Fetches the PDB content for a specific result.
    """

    def __init__(self, result: SearchResult):
        super().__init__()
        self.result = result
        self.signals = WorkerSignals()

    @Slot()
    def run(self):
        try:
            # 1. Cache
            if self.result.source == "Cache":
                pdb = MoleculeService.get_cached_structure(self.result.identifier)
                if pdb:
                    self.signals.finished.emit(pdb)
                    return
                # If source was cache but file missing, fallback (shouldn't happen logic-wise but good for safety)

            # 2. PubChem
            if self.result.source == "PubChem":
                # Convert logic
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{self.result.identifier}/SDF"
                res = requests.get(url, timeout=10)
                if res.status_code == 200:
                    pdb = MoleculeService.convert_smiles_to_pdb(res.text)
                    save_name = self.result.inquiry_name if self.result.inquiry_name else self.result.identifier
                    MoleculeService.save_to_cache(save_name, pdb)
                    self.signals.finished.emit(pdb)
                    return

            # 3. PDB
            if self.result.source == "PDB":
                url = f"https://files.rcsb.org/download/{self.result.identifier}.pdb"
                res = requests.get(url, timeout=10)
                if res.status_code == 200:
                    pdb = res.text
                    save_name = self.result.inquiry_name if self.result.inquiry_name else self.result.identifier
                    MoleculeService.save_to_cache(save_name, pdb)
                    self.signals.finished.emit(pdb)
                    return

            raise ValueError(f"Could not fetch data for {self.result.identifier} from {self.result.source}")

        except Exception as e:
            self.signals.error.emit(str(e))
