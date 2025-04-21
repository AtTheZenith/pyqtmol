# pyqtmol

[![Python ≥3.12](https://img.shields.io/badge/python-%3E%3D3.12-blue)](https://www.python.org/downloads/)  
[![License: GPL‑3.0](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

An **interactive 3D compound structure visualizer** built in Python.  
Under the hood it uses:
- **PySide6** (Qt) for a modern desktop GUI  
- **3Dmol.js** for high‑quality WebGL rendering  
- **RDKit** to generate 3D conformers from SMILES  
- **PubChem** & **RCSB PDB** as structure data sources  

---

## 🌟 Key Features

- **Search by name**  
  Fetches SMILES from PubChem and falls back to PDB from RCSB if needed.  
- **Native WebView**  
  Renders within a Qt `QWebEngineView` using the 3Dmol.js API.  
- **Dynamic caching**  
  Caches downloaded `.pdb` files and the 3Dmol.js bundle for low data usage.  

---

## 🚀 Installation

### Using `pip`

1.  **Clone and enter folder**  
    ```bash
    git clone https://github.com/AtTheZenith/pyqtmol.git
    cd pyqtmol
    ```

2.  **(Optional) Create a virtual environment**
    ```bash
    python -m venv
    # Linux/macOS
    source .venv/bin/activate
    # Windows
    .venv\Scripts\activate
    ```

3.  **Install dependencies**
    ```bash
    pip install --upgrade pip
    pip install -r requirements.txt
    ```

4.  **Run the app**
    ```bash
    python main.py
    ```

### Using `uv`

1. **Clone and enter**
    ```bash
    git clone https://github.com/AtTheZenith/pyqtmol.git
    cd pyqtmol
    ```

2. **Create and enter a virtual enviroment**
    ```bash
    uv venv
    .venv/Scripts/activate
    ```

3. **Run the app**
    ```bash
    uv run main.py
    ```
