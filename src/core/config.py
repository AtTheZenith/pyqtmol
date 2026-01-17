import os
import sys

def get_base_path() -> str:
    """Returns the base path of the application."""
    # Handle the case where we are running from a frozen bundle (e.g. PyInstaller)
    if getattr(sys, 'frozen', False):
         return sys._MEIPASS # type: ignore
    return os.getcwd()

CACHE_DIR = os.path.join(get_base_path(), "cache")
ASSETS_DIR = os.path.join(get_base_path(), "assets")

def ensure_dirs():
    """Ensure critical directories exist."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    os.makedirs(ASSETS_DIR, exist_ok=True)
