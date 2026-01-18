import os

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# Assets are now inside src/resources
# config.py is in src/core, so resources is ../resources
SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSETS_DIR = os.path.join(SRC_DIR, "resources")
TEMPLATES_DIR = os.path.join(SRC_DIR, "templates")
