import sys
from PySide6.QtWidgets import QApplication

from src.core.config import ensure_dirs
from src.core.assets import check_and_download_assets
from src.ui.window import MainWindow

def main():
    # 1. Setup Environment
    ensure_dirs()
    
    # 2. Check Assets (Blocking here is acceptable as it's startup, 
    # but ideally show a splash screen. For now, CLI log is fine or fast enough)
    check_and_download_assets()

    # 3. Launch App
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
