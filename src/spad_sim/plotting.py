import numpy as np
import matplotlib.pyplot as plt

def plot_coincidence_histogram(time_differences_s, window_ns=5.0, coincidence_window_ns=0.1, output_file=None):
    """
    Calculates and plots a histogram of time differences.

    The number of bins is now calculated based on the desired coincidence window.

    Args:
        time_differences_s (np.ndarray): Array of time differences in seconds.
        window_ns (float): The total time to show on the histogram, from -window_ns to +window_ns.
        coincidence_window_ns (float): The desired width of each bin in nanoseconds.
        output_file (str, optional): If provided, the plot will be saved to this file. Defaults to None.
    """
    print("\nPlotting coincidence histogram...")

    # --- NEW: Calculate the number of bins ---
    # The total span of the histogram is 2 * window_ns
    total_histogram_span_ns = 2 * window_ns
    num_bins = int(total_histogram_span_ns / coincidence_window_ns)

    print(f"Histogram Bin Width (Coincidence Window): {coincidence_window_ns} ns")
    print(f"Calculated Number of Bins: {num_bins}")

    # Convert differences to nanoseconds for plotting
    time_differences_ns = time_differences_s * 1e9

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(time_differences_ns, bins=num_bins, range=(-window_ns, window_ns), color='seagreen')
    plt.title('Coincidence Plot', fontsize=16)
    plt.xlabel('Time Difference (τ = t₂ - t₁) [ns]', fontsize=12)
    plt.ylabel('Coincidence Counts', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)

    if output_file:
        print(f"Saving plot to {output_file}...")
        plt.savefig(output_file)
        plt.close()
    else:
        plt.show()
