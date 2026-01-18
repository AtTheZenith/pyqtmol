from PySide6.QtWidgets import QListWidget, QStyledItemDelegate, QListWidgetItem, QStyle
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QPainter, QColor, QFont

from src.ui.theme import Theme


class ResultDelegate(QStyledItemDelegate):
    def paint(self, painter: QPainter, option, index):
        painter.save()
        painter.setRenderHint(QPainter.Antialiasing)

        rect = option.rect
        data = index.data(Qt.UserRole)

        # Background
        if option.state & QStyle.State_Selected:
            painter.translate(2, 0)
            bg_color = QColor(Theme.SURFACE_HOVER)
            # Active Indicator
            painter.setPen(Qt.NoPen)
            painter.setBrush(QColor(Theme.PRIMARY))
            painter.drawRoundedRect(rect.left(), rect.top() + 10, 4, rect.height() - 20, 2, 2)
        elif option.state & QStyle.State_MouseOver:
            bg_color = QColor(Theme.SURFACE)
        else:
            bg_color = Qt.transparent

        painter.setPen(Qt.NoPen)
        painter.setBrush(bg_color)
        painter.drawRoundedRect(rect.adjusted(8, 2, -8, -2), 8, 8)

        # Text
        painter.setPen(QColor(Theme.TEXT_PRIMARY))
        title_font = QFont("Segoe UI", 11, QFont.Bold)
        painter.setFont(title_font)

        title_rect = rect.adjusted(20, 10, -10, -25)
        painter.drawText(title_rect, Qt.AlignLeft | Qt.AlignTop, data.name)

        # Subtext
        painter.setPen(QColor(Theme.TEXT_SECONDARY))
        sub_font = QFont("Segoe UI", 9)
        painter.setFont(sub_font)

        sub_rect = rect.adjusted(20, 32, -10, -10)
        painter.drawText(sub_rect, Qt.AlignLeft | Qt.AlignTop, f"Source: {data.source}")

        painter.restore()

    def sizeHint(self, option, index):
        return QSize(option.rect.width(), 65)


class ResultList(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameShape(QListWidget.NoFrame)
        self.setObjectName("result_list")
        self.setItemDelegate(ResultDelegate())
        self.setVerticalScrollMode(QListWidget.ScrollPerPixel)

    def add_result(self, result):
        item = QListWidgetItem()
        item.setData(Qt.UserRole, result)
        self.addItem(item)
