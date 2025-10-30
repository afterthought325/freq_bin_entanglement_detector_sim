import sys
import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import QThread, pyqtSignal
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QPushButton, QGroupBox, QFormLayout, QLabel, QLineEdit, QProgressBar
)
from .core import simulate_spad_time_tags_iterative, find_time_differences

class SimulationWorker(QThread):
    progress = pyqtSignal(int)
    new_data = pyqtSignal(object)
    finished = pyqtSignal()

    def __init__(self, params, search_window_s):
        super().__init__()
        self.params = params
        self.search_window_s = search_window_s
        self.running = True

    def run(self):
        num_chunks = 10
        all_tags1, all_tags2 = [], []

        for i, (tags1_chunk, tags2_chunk) in enumerate(simulate_spad_time_tags_iterative(num_chunks=num_chunks, **self.params)):
            if not self.running:
                break

            all_tags1.extend(tags1_chunk)
            all_tags2.extend(tags2_chunk)

            time_diffs = find_time_differences(np.array(all_tags1), np.array(all_tags2), self.search_window_s)
            self.new_data.emit(time_diffs)
            self.progress.emit(int((i + 1) / num_chunks * 100))

        self.finished.emit()

    def stop(self):
        self.running = False


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("SPAD Simulation Tool")

        # Main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)

        # Left panel for inputs and controls
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        main_layout.addWidget(left_panel)

        self._create_input_parameters_group(left_layout)

        self.start_button = QPushButton("Start")
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.setEnabled(False)

        self.start_button.clicked.connect(self.start_simulation)
        self.cancel_button.clicked.connect(self.cancel_simulation)

        left_layout.addWidget(self.start_button)
        left_layout.addWidget(self.cancel_button)
        self.progress_bar = QProgressBar()
        left_layout.addWidget(self.progress_bar)
        left_layout.addStretch()


        # Right panel for outputs
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)
        main_layout.addWidget(right_panel)

        self._create_output_display(right_layout)

    def _create_output_display(self, parent_layout):
        # Histogram Plot
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')
        self.plot_widget.setTitle("Coincidence Histogram")
        self.plot_widget.setLabel('left', 'Counts')
        self.plot_widget.setLabel('bottom', 'Time Difference (ns)')
        parent_layout.addWidget(self.plot_widget)

        # Stats Panel
        stats_group = QGroupBox("Simulation Stats")
        stats_layout = QFormLayout(stats_group)
        self.stats = {
            'peak_offset': QLabel("N/A"),
            'fwhm': QLabel("N/A"),
            'std_dev': QLabel("N/A"),
        }
        stats_layout.addRow(QLabel("Peak Offset (ns):"), self.stats['peak_offset'])
        stats_layout.addRow(QLabel("FWHM (ns):"), self.stats['fwhm'])
        stats_layout.addRow(QLabel("Std Dev (ns):"), self.stats['std_dev'])
        parent_layout.addWidget(stats_group)

    def start_simulation(self):
        self.start_button.setEnabled(False)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setValue(0)
        self.plot_widget.clear()
        for stat in self.stats.values():
            stat.setText("N/A")

        sim_params = {key: float(widget.text()) for key, widget in self.params.items()}

        # Using a fixed value for now, can be made configurable
        search_window_s = 0.1

        self.worker = SimulationWorker(sim_params, search_window_s)
        self.worker.progress.connect(self.update_progress)
        self.worker.new_data.connect(self.update_plot)
        self.worker.finished.connect(self.simulation_finished)
        self.worker.start()

    def cancel_simulation(self):
        if self.worker:
            self.worker.stop()
            self.cancel_button.setEnabled(False)
            self.start_button.setEnabled(True)

    def update_progress(self, value):
        self.progress_bar.setValue(value)

    def update_plot(self, time_diffs):
        if len(time_diffs) > 1:
            bins = np.linspace(-5, 5, 200) # Example binning
            y, x = np.histogram(time_diffs * 1e9, bins=bins)
            self.plot_widget.plot(x, y, stepMode="center", fillLevel=0, brush=(0,0,255,150), clear=True)
            self._calculate_and_display_stats(time_diffs)

    def simulation_finished(self):
        self.start_button.setEnabled(True)
        self.cancel_button.setEnabled(False)
        self.progress_bar.setValue(100)

    def _calculate_and_display_stats(self, time_diffs):
        if len(time_diffs) < 2:
            return

        # Peak Offset
        hist, bins = np.histogram(time_diffs, bins=100)
        peak_bin = np.argmax(hist)
        peak_offset_s = (bins[peak_bin] + bins[peak_bin+1]) / 2
        self.stats['peak_offset'].setText(f"{peak_offset_s * 1e9:.3f}")

        # Standard Deviation
        std_dev_s = np.std(time_diffs)
        self.stats['std_dev'].setText(f"{std_dev_s * 1e9:.3f}")

        # FWHM (Full Width at Half Maximum)
        try:
            half_max = np.max(hist) / 2.0
            above_half_max = np.where(hist > half_max)[0]
            first_idx, last_idx = above_half_max[0], above_half_max[-1]
            fwhm = bins[last_idx] - bins[first_idx]
            self.stats['fwhm'].setText(f"{fwhm * 1e9:.3f}")
        except IndexError:
            self.stats['fwhm'].setText("N/A")

    def _create_input_parameters_group(self, parent_layout):
        group_box = QGroupBox("Input Parameters")
        form_layout = QFormLayout(group_box)

        # Default parameters from spad_sim.py
        self.params = {
            'simulation_time_s': QLineEdit("1.0"),
            'source_pair_rate_hz': QLineEdit("5e5"),
            'distance_km': QLineEdit("5.0"),
            'dispersion_ps_nm_km': QLineEdit("17.0"),
            'path_loss1_db': QLineEdit("10.0"),
            'path_loss2_db': QLineEdit("10.0"),
            'wavelength1_nm': QLineEdit("1530"),
            'qe1': QLineEdit("0.25"),
            'dcr1_hz': QLineEdit("600000"),
            'jitter1_s': QLineEdit("200e-12"),
            'dead_time1_s': QLineEdit("3e-6"),
            'wavelength2_nm': QLineEdit("1550"),
            'qe2': QLineEdit("0.25"),
            'dcr2_hz': QLineEdit("600000"),
            'jitter2_s': QLineEdit("150e-12"),
            'dead_time2_s': QLineEdit("3e-6"),
        }

        # Simulation Parameters
        sim_group = QGroupBox("Simulation")
        sim_layout = QFormLayout(sim_group)
        sim_layout.addRow(QLabel("Simulation Time (s):"), self.params['simulation_time_s'])
        sim_layout.addRow(QLabel("Source Pair Rate (Hz):"), self.params['source_pair_rate_hz'])
        sim_layout.addRow(QLabel("Distance (km):"), self.params['distance_km'])
        sim_layout.addRow(QLabel("Dispersion (ps/nm/km):"), self.params['dispersion_ps_nm_km'])
        sim_layout.addRow(QLabel("Path Loss 1 (dB):"), self.params['path_loss1_db'])
        sim_layout.addRow(QLabel("Path Loss 2 (dB):"), self.params['path_loss2_db'])
        form_layout.addWidget(sim_group)

        # Detector 1 Parameters
        det1_group = QGroupBox("Detector 1")
        det1_layout = QFormLayout(det1_group)
        det1_layout.addRow(QLabel("Wavelength (nm):"), self.params['wavelength1_nm'])
        det1_layout.addRow(QLabel("Quantum Efficiency:"), self.params['qe1'])
        det1_layout.addRow(QLabel("Dark Count Rate (Hz):"), self.params['dcr1_hz'])
        det1_layout.addRow(QLabel("Timing Jitter (s):"), self.params['jitter1_s'])
        det1_layout.addRow(QLabel("Dead Time (s):"), self.params['dead_time1_s'])
        form_layout.addWidget(det1_group)

        # Detector 2 Parameters
        det2_group = QGroupBox("Detector 2")
        det2_layout = QFormLayout(det2_group)
        det2_layout.addRow(QLabel("Wavelength (nm):"), self.params['wavelength2_nm'])
        det2_layout.addRow(QLabel("Quantum Efficiency:"), self.params['qe2'])
        det2_layout.addRow(QLabel("Dark Count Rate (Hz):"), self.params['dcr2_hz'])
        det2_layout.addRow(QLabel("Timing Jitter (s):"), self.params['jitter2_s'])
        det2_layout.addRow(QLabel("Dead Time (s):"), self.params['dead_time2_s'])
        form_layout.addWidget(det2_group)

        parent_layout.addWidget(group_box)

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
