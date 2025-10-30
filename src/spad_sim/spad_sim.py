import argparse
from .core import simulate_spad_time_tags, find_time_differences
from .io import save_tags_to_csv
from .plotting import plot_coincidence_histogram

def main():
    parser = argparse.ArgumentParser(description="SPAD Simulation Tool")

    # Simulation parameters
    parser.add_argument('--simulation_time_s', type=float, default=1.0, help='Simulation time in seconds')
    parser.add_argument('--source_pair_rate_hz', type=float, default=5e5, help='Source pair rate in Hz')
    parser.add_argument('--distance_km', type=float, default=5.0, help='Distance in kilometers')
    parser.add_argument('--dispersion_ps_nm_km', type=float, default=17.0, help='Dispersion in ps/nm/km')
    parser.add_argument('--path_loss1_db', type=float, default=10.0, help='Path loss 1 in dB')
    parser.add_argument('--path_loss2_db', type=float, default=10.0, help='Path loss 2 in dB')

    # Detector 1 parameters
    parser.add_argument('--wavelength1_nm', type=float, default=1530, help='Wavelength for detector 1 in nm')
    parser.add_argument('--qe1', type=float, default=0.25, help='Quantum efficiency for detector 1')
    parser.add_argument('--dcr1_hz', type=float, default=600000, help='Dark count rate for detector 1 in Hz')
    parser.add_argument('--jitter1_s', type=float, default=200e-12, help='Timing jitter for detector 1 in seconds')
    parser.add_argument('--dead_time1_s', type=float, default=3e-6, help='Dead time for detector 1 in seconds')

    # Detector 2 parameters
    parser.add_argument('--wavelength2_nm', type=float, default=1550, help='Wavelength for detector 2 in nm')
    parser.add_argument('--qe2', type=float, default=0.25, help='Quantum efficiency for detector 2')
    parser.add_argument('--dcr2_hz', type=float, default=600000, help='Dark count rate for detector 2 in Hz')
    parser.add_argument('--jitter2_s', type=float, default=150e-12, help='Timing jitter for detector 2 in seconds')
    parser.add_argument('--dead_time2_s', type=float, default=3e-6, help='Dead time for detector 2 in seconds')

    # Post-processing and plotting parameters
    parser.add_argument('--output_file', type=str, default="time_tags.csv", help='Output CSV file name')
    parser.add_argument('--no-save', action='store_true', help="Don't save the time tags to a file")
    parser.add_argument('--plot-window-ns', type=float, default=10.0, help='Plotting window in nanoseconds')
    parser.add_argument('--coincidence-window-ns', type=float, default=0.08, help='Coincidence window for histogram bins in nanoseconds')
    parser.add_argument('--no-plot', action='store_true', help="Don't show the coincidence plot")
    parser.add_argument('--gui', action='store_true', help="Launch the GUI")


    args = parser.parse_args()

    if args.gui:
        from . import gui
        gui.main()
        return

    # Run the simulation
    sim_params = {
        'simulation_time_s': args.simulation_time_s,
        'source_pair_rate_hz': args.source_pair_rate_hz,
        'distance_km': args.distance_km,
        'dispersion_ps_nm_km': args.dispersion_ps_nm_km,
        'path_loss1_db': args.path_loss1_db,
        'path_loss2_db': args.path_loss2_db,
        'wavelength1_nm': args.wavelength1_nm,
        'qe1': args.qe1,
        'dcr1_hz': args.dcr1_hz,
        'jitter1_s': args.jitter1_s,
        'dead_time1_s': args.dead_time1_s,
        'wavelength2_nm': args.wavelength2_nm,
        'qe2': args.qe2,
        'dcr2_hz': args.dcr2_hz,
        'jitter2_s': args.jitter2_s,
        'dead_time2_s': args.dead_time2_s,
    }
    time_tags_1, time_tags_2 = simulate_spad_time_tags(**sim_params)

    # Save the results
    if not args.no_save:
        save_tags_to_csv(time_tags_1, time_tags_2, filename=args.output_file)

    # Plot the results
    if not args.no_plot:
        search_window_s = args.plot_window_ns * 1e-9
        time_differences = find_time_differences(time_tags_1, time_tags_2, search_window_s)
        plot_coincidence_histogram(
            time_differences,
            window_ns=args.plot_window_ns,
            coincidence_window_ns=args.coincidence_window_ns
        )

if __name__ == '__main__':
    main()
