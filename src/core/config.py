import os

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)

# config.py is in src/core, so resources is ../resources
SRC_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESOURCES_DIR = os.path.join(SRC_DIR, "resources")
TEMPLATES_DIR = os.path.join(SRC_DIR, "templates")

if not os.path.exists(RESOURCES_DIR):
    RESOURCES_DIR = CACHE_DIR

if not os.path.exists(TEMPLATES_DIR):
    TEMPLATES_DIR = CACHE_DIR
