import os
from rdkit import Chem
from rdkit.Chem import rdDistGeom

from src.core.config import CACHE_DIR


class MoleculeService:
    """
    Service for handling molecule data conversions and file operations.
    Keeps IO and heavy logic separate from UI.
    """

    @staticmethod
    def get_cached_structure(compound_name: str) -> str | None:
        """Retrieves the PDB structure from cache if it exists."""
        safe_name = compound_name.upper().strip()
        path = os.path.join(CACHE_DIR, f"{safe_name}.pdb")
        if os.path.isfile(path):
            with open(path, "r") as f:
                return f.read()
        return None

    @staticmethod
    def save_to_cache(compound_name: str, pdb_data: str) -> str:
        """Saves the PDB structure to cache."""
        safe_name = compound_name.upper().strip()
        path = os.path.join(CACHE_DIR, f"{safe_name}.pdb")
        with open(path, "w") as f:
            f.write(pdb_data)
        return path

    @staticmethod
    def convert_smiles_to_pdb(mol_block: str) -> str:
        """
        Converts a MolBlock (from SDF/SMILES) to a PDB block with 3D coordinates.
        This is a CPU-intensive operation.
        """
        mol = Chem.MolFromMolBlock(mol_block, sanitize=False, removeHs=False)
        if mol is None:
            raise ValueError("Could not parse molecule block")

        Chem.SanitizeMol(mol)
        rdDistGeom.EmbedMolecule(mol, rdDistGeom.ETKDG())
        return Chem.MolToPDBBlock(mol)
