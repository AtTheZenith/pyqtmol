import os
import requests

from src.core.config import CACHE_DIR, RESOURCES_DIR, TEMPLATES_DIR


REQUIRED_ASSETS = {
    os.path.join(CACHE_DIR, "3Dmol-min.js"): "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.1/3Dmol-min.js",
    os.path.join(RESOURCES_DIR, "check.svg"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/src/resources/check.svg",
    os.path.join(RESOURCES_DIR, "icon.ico"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/src/resources/icon.ico",
    os.path.join(RESOURCES_DIR, "style.qss"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/src/resources/style.qss",
    os.path.join(TEMPLATES_DIR, "none.html"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/src/templates/none.html",
    os.path.join(TEMPLATES_DIR, "template.html"): "https://raw.githubusercontent.com/AtTheZenith/pyqtmol/main/src/templates/template.html",
}


def integrity_check():
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
