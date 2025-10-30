import sys
import os
from pyvirtualdisplay import Display
from PyQt6.QtWidgets import QApplication
from PyQt6.QtCore import QTimer
from spad_sim.gui import MainWindow

def main():
    with Display():
        app = QApplication(sys.argv)
        window = MainWindow()
        window.show()

        # Take a screenshot after a short delay to allow the window to render
        QTimer.singleShot(1000, lambda: take_screenshot(window))

        # Exit after a slightly longer delay
        QTimer.singleShot(1500, app.quit)

        app.exec()

def take_screenshot(window):
    screenshot = window.grab()
    screenshot.save("jules-scratch/verification/verification.png", "png")
    print("Screenshot saved to jules-scratch/verification/verification.png")

if __name__ == '__main__':
    # Add src to the python path to allow the import of spad_sim
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'src')))
    main()
