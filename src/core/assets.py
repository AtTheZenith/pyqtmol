import os
import requests
from src.core.config import ASSETS_DIR, CACHE_DIR

REQUIRED_ASSETS = {
    os.path.join(ASSETS_DIR, "64x.ico"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/64x.ico",
    os.path.join(ASSETS_DIR, "none.html"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/none.html",
    os.path.join(ASSETS_DIR, "styles.qss"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/styles.qss",
    os.path.join(ASSETS_DIR, "template.html"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/assets/template.html",
    os.path.join(CACHE_DIR, "3Dmol-min.js"): "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.1/3Dmol-min.js",
}

def check_and_download_assets():
    """Checks for required assets and downloads them if missing."""
    print("Checking assets...")
    for file_path, url in REQUIRED_ASSETS.items():
        if not os.path.exists(file_path):
            print(f"Downloading {os.path.basename(file_path)}...")
            try:
                res = requests.get(url, timeout=30)
                if res.status_code == 200:
                    with open(file_path, "wb") as f:
                        f.write(res.content)
                else:
                    print(f"Failed to download {url}: Status {res.status_code}")
            except Exception as e:
                print(f"Error downloading {url}: {e}")
