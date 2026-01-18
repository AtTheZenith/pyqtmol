import sys

from PySide6.QtWidgets import QApplication

from src.core.assets import integrity_check

from src.ui.window import MainWindow


def main():
    integrity_check()

    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
