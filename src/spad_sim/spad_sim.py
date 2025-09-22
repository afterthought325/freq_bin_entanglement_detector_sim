import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def simulate_spad_time_tags(
    # New parameters for loss
    source_pair_rate_hz=5e6,
    path_loss1_db=0.0,
    path_loss2_db=0.0,
    # Existing parameters
    distance_km=0.0,
    wavelength1_nm=1530,
    wavelength2_nm=1550,
    dispersion_ps_nm_km=17.0,
    simulation_time_s=0.1,
    qe1=0.7,
    dcr1_hz=100,
    jitter1_s=50e-12,
    dead_time1_s=50e-9,
    qe2=0.75,
    dcr2_hz=120,
    jitter2_s=60e-12,
    dead_time2_s=45e-9,
):
    """
    Simulates time tags, now including optical path loss.
    """
    print("Starting simulation...")

    # --- 1. Calculate Transmission Probability from Path Loss ---
    # This converts the dB loss into a linear survival probability for each photon.
    transmission1 = 10**(-path_loss1_db / 10)
    transmission2 = 10**(-path_loss2_db / 10)
    print(f"Path 1 Loss: {path_loss1_db} dB -> {transmission1:.3f} Transmission")
    print(f"Path 2 Loss: {path_loss2_db} dB -> {transmission2:.3f} Transmission")

    # --- 2. Generate Photon Pair Creation Times ---
    num_pairs = int(simulation_time_s * source_pair_rate_hz)
    pair_intervals = np.random.exponential(1.0 / source_pair_rate_hz, num_pairs)
    pair_creation_times = np.cumsum(pair_intervals)
    pair_creation_times = pair_creation_times[pair_creation_times < simulation_time_s]
    
    print(f"Generated {len(pair_creation_times)} photon pairs at the source.")

    # --- 3. Calculate Chromatic Dispersion Delay ---
    wavelength_diff_nm = wavelength2_nm - wavelength1_nm
    delay_ps = dispersion_ps_nm_km * distance_km * wavelength_diff_nm
    delay_s = delay_ps * 1e-12
    
    if distance_km > 0:
        print(f"Calculated Dispersion Delay: {delay_ps:.2f} ps ({delay_s * 1e9:.2f} ns)")

    # --- 4. Simulate Photon Detections ---
    # The detection probability is now the path transmission * the detector QE.
    # This single check efficiently models both effects.
    
    # Detector 1
    total_efficiency1 = transmission1 * qe1
    detection_mask1 = np.random.rand(len(pair_creation_times)) < total_efficiency1
    photon_times1 = pair_creation_times[detection_mask1]

    # Detector 2
    total_efficiency2 = transmission2 * qe2
    detection_mask2 = np.random.rand(len(pair_creation_times)) < total_efficiency2
    photon_times2 = pair_creation_times[detection_mask2]

    # --- 5. Generate Dark Count Times ---
    num_dark_counts1 = np.random.poisson(dcr1_hz * simulation_time_s)
    dark_count_times1 = np.random.uniform(0, simulation_time_s, num_dark_counts1)

    num_dark_counts2 = np.random.poisson(dcr2_hz * simulation_time_s)
    dark_count_times2 = np.random.uniform(0, simulation_time_s, num_dark_counts2)

    # --- 6. Combine, Delay, and Add Jitter ---
    all_events1 = np.concatenate((photon_times1, dark_count_times1))
    all_events1 += np.random.normal(0, jitter1_s, len(all_events1))
    all_events1.sort()

    all_events2 = np.concatenate((photon_times2, dark_count_times2))
    all_events2 += delay_s  # Apply dispersion delay
    all_events2 += np.random.normal(0, jitter2_s, len(all_events2))
    all_events2.sort()

    # --- 7. Apply Detector Dead Time ---
    def apply_dead_time(events, dead_time):
        if len(events) == 0: return np.array([])
        final_tags, last_detection_time = [], -np.inf
        for t in events:
            if (t - last_detection_time) > dead_time:
                final_tags.append(t)
                last_detection_time = t
        return np.array(final_tags)

    print("Applying dead time...")
    final_tags1 = apply_dead_time(all_events1, dead_time1_s)
    final_tags2 = apply_dead_time(all_events2, dead_time2_s)
    
    print("Simulation finished!")
    print(f"Detector 1 registered {len(final_tags1)} events.")
    print(f"Detector 2 registered {len(final_tags2)} events.")

    return final_tags1, final_tags2

def save_tags_to_csv(tags1, tags2, filename="time_tags.csv"):
    """Saves the two time tag arrays to a single CSV file."""
    print(f"\nSaving time tags to {filename}...")
    # Use pandas to easily handle the two arrays of different lengths
    df = pd.DataFrame({
        'Detector1_Time_s': pd.Series(tags1),
        'Detector2_Time_s': pd.Series(tags2)
    })
    df.to_csv(filename, index=False)
    print("Save complete.")

# --- NEW Coincidence Counting Function ---
def find_time_differences(tags1, tags2, search_window_s):
    """
    Finds all time differences between two tag streams within a search window.
    This is the core data processing function.

    Args:
        tags1 (np.ndarray): Time tags from detector 1.
        tags2 (np.ndarray): Time tags from detector 2.
        search_window_s (float): The +/- window in seconds to look for a match.

    Returns:
        np.ndarray: An array of all found time differences (t2 - t1) in seconds.
    """
    print(f"\nSearching for coincidences within a +/- {search_window_s * 1e9:.2f} ns window...")
    time_differences = []
    idx2 = 0

    for t1 in tags1:
        # Advance the second pointer to the start of the search window
        while idx2 < len(tags2) and tags2[idx2] < t1 - search_window_s:
            idx2 += 1

        # Check all tags that fall within the window [t1 - W, t1 + W]
        temp_idx2 = idx2
        while temp_idx2 < len(tags2) and tags2[temp_idx2] <= t1 + search_window_s:
            diff = tags2[temp_idx2] - t1
            time_differences.append(diff)
            temp_idx2 += 1

    print(f"Found {len(time_differences)} potential coincidence events.")
    return np.array(time_differences)

def plot_coincidence_histogram(tags1, tags2, window_ns=5.0, coincidence_window_ns=0.1):
    """
    Calculates and plots a histogram of time differences.

    The number of bins is now calculated based on the desired coincidence window.

    Args:
        tags1 (np.ndarray): Time tags from detector 1.
        tags2 (np.ndarray): Time tags from detector 2.
        window_ns (float): The total time to show on the histogram, from -window_ns to +window_ns.
        coincidence_window_ns (float): The desired width of each bin in nanoseconds.
    """
    print("\nCalculating coincidences...")
    time_differences = []

    # This search window must be the same as the histogram window to find all relevant events
    search_window_s = window_ns * 1e-9
    idx2 = 0

    for t1 in tags1:
        while idx2 < len(tags2) and tags2[idx2] < t1 - search_window_s:
            idx2 += 1
        temp_idx2 = idx2
        while temp_idx2 < len(tags2) and tags2[temp_idx2] <= t1 + search_window_s:
            time_differences.append(tags2[temp_idx2] - t1)
            temp_idx2 += 1

    print(f"Found {len(time_differences)} raw coincidences within a +/- {window_ns} ns window.")

    # --- NEW: Calculate the number of bins ---
    # The total span of the histogram is 2 * window_ns
    total_histogram_span_ns = 2 * window_ns
    num_bins = int(total_histogram_span_ns / coincidence_window_ns)

    print(f"Histogram Bin Width (Coincidence Window): {coincidence_window_ns} ns")
    print(f"Calculated Number of Bins: {num_bins}")

    # Convert differences to nanoseconds for plotting
    time_differences_ns = np.array(time_differences) * 1e9

    # Plotting the histogram
    plt.figure(figsize=(10, 6))
    plt.hist(time_differences_ns, bins=num_bins, range=(-window_ns, window_ns), color='seagreen')
    plt.title('Coincidence Plot', fontsize=16)
    plt.xlabel('Time Difference (τ = t₂ - t₁) [ns]', fontsize=12)
    plt.ylabel('Coincidence Counts', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.show()


if __name__ == '__main__':
    # --- How to use the simulation ---
    
    # You can adjust these parameters to match your experimental setup
    params = {
        'simulation_time_s': 1.0,       # Simulate for 100 milliseconds
        'source_pair_rate_hz': 2.5e6, # Let's use a very bright 50 MHz source
        'distance_km': 05.0,             # 10 kilometers of fiber
        'dispersion_ps_nm_km': 17.0,     # Standard SMF-28 fiber dispersion
        'path_loss1_db': 12.0,       # 10 dB loss (90% lost) in path 1
        'path_loss2_db': 12.0,       # 13 dB loss (~95% lost) in path 2
        'wavelength1_nm': 1530,          # Wavelength for detector 1
        'qe1': 0.2,                    # 65% QE for detector 1
        'dcr1_hz': 600000,                 # 150 dark counts/sec for detector 1
        'jitter1_s': 200e-12,            # 40 picosecond jitter for detector 1
        'dead_time1_s': 3e-6,          # 60 nanosecond dead time for detector 1
        'wavelength2_nm': 1550,          # Wavelength for detector 2
        'qe2': 0.2,                    # 70% QE for detector 2
        'dcr2_hz': 600000,                 # 200 dark counts/sec for detector 2
        'jitter2_s': 150e-12,            # 45 picosecond jitter for detector 2
        'dead_time2_s': 3e-6,          # 55 nanosecond dead time for detector 2
    }

    # Run the simulation
    time_tags_1, time_tags_2 = simulate_spad_time_tags(**params)

    # 1. Save the results to a CSV file
    save_tags_to_csv(time_tags_1, time_tags_2)

    # Now we specify the total window and the desired resolution (bin width).
    # A bin width of 0.08 ns (80 ps) is a good choice, as it's around twice the jitter.
    plot_coincidence_histogram(
        time_tags_1,
        time_tags_2,
        window_ns=20.0,
        coincidence_window_ns=0.16
    )
