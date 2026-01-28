import os

from PySide6.QtGui import QColor, QPalette
from PySide6.QtWidgets import QApplication


from src.core.config import RESOURCES_DIR


class Theme:
    # Colors

    BACKGROUND = "#1e1e1e"
    SURFACE = "#252526"
    SURFACE_HOVER = "#2a2d2e"
    PRIMARY = "#3b82f6"
    PRIMARY_HOVER = "#60a5fa"
    TEXT_PRIMARY = "#ffffff"
    TEXT_SECONDARY = "#b0b0b0"
    BORDER = "#3e3e42"
    ERROR = "#ef4444"
    SUCCESS = "#22c55e"

    @staticmethod
    def apply(app: QApplication):
        app.setStyle("Fusion")
        palette = QPalette()
        palette.setColor(QPalette.Window, QColor(Theme.BACKGROUND))
        palette.setColor(QPalette.WindowText, QColor(Theme.TEXT_PRIMARY))
        palette.setColor(QPalette.Base, QColor(Theme.SURFACE))
        palette.setColor(QPalette.AlternateBase, QColor(Theme.SURFACE_HOVER))
        palette.setColor(QPalette.ToolTipBase, QColor(Theme.SURFACE))
        palette.setColor(QPalette.ToolTipText, QColor(Theme.TEXT_PRIMARY))
        palette.setColor(QPalette.Text, QColor(Theme.TEXT_PRIMARY))
        palette.setColor(QPalette.Button, QColor(Theme.SURFACE))
        palette.setColor(QPalette.ButtonText, QColor(Theme.TEXT_PRIMARY))
        palette.setColor(QPalette.BrightText, QColor(Theme.ERROR))
        palette.setColor(QPalette.Link, QColor(Theme.PRIMARY))
        palette.setColor(QPalette.Highlight, QColor(Theme.PRIMARY))
        palette.setColor(QPalette.HighlightedText, QColor(Theme.TEXT_PRIMARY))
        app.setPalette(palette)

        # Load QSS

        qss_path = os.path.join(RESOURCES_DIR, "style.qss")

        if os.path.exists(qss_path):
            with open(qss_path, "r") as f:
                qss = f.read()

            # Replace placeholders
            # We iterate over class attributes to find replacements
            vars = {
                "@BACKGROUND": Theme.BACKGROUND,
                "@SURFACE": Theme.SURFACE,
                "@SURFACE_HOVER": Theme.SURFACE_HOVER,
                "@PRIMARY": Theme.PRIMARY,
                "@PRIMARY_HOVER": Theme.PRIMARY_HOVER,
                "@TEXT_PRIMARY": Theme.TEXT_PRIMARY,
                "@TEXT_SECONDARY": Theme.TEXT_SECONDARY,
                "@BORDER": Theme.BORDER,
                "@ERROR": Theme.ERROR,
                "@SUCCESS": Theme.SUCCESS,
                "@CHECK": os.path.join(RESOURCES_DIR, "check.svg").replace("\\", "/"),
            }

            for key in sorted(vars.keys(), key=len, reverse=True):
                val = vars[key]
                qss = qss.replace(key, val)

            app.setStyleSheet(qss)
        else:
            print(f"Warning: Stylesheet not found at {qss_path}")
