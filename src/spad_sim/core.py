import numpy as np

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


def _apply_dead_time_iterative(events, dead_time, last_detection_time):
    """Helper function to apply dead time while preserving state."""
    if len(events) == 0:
        return np.array([]), last_detection_time

    final_tags = []
    for t in events:
        if (t - last_detection_time) > dead_time:
            final_tags.append(t)
            last_detection_time = t
    return np.array(final_tags), last_detection_time


def simulate_spad_time_tags_iterative(
    num_chunks=10,
    # Pass all other simulation parameters
    **kwargs
):
    """
    A generator that simulates time tags in chunks, yielding intermediate results.
    """
    simulation_time_s = kwargs.get('simulation_time_s', 0.1)
    source_pair_rate_hz = kwargs.get('source_pair_rate_hz', 5e6)
    path_loss1_db = kwargs.get('path_loss1_db', 0.0)
    path_loss2_db = kwargs.get('path_loss2_db', 0.0)
    distance_km = kwargs.get('distance_km', 0.0)
    wavelength1_nm = kwargs.get('wavelength1_nm', 1530)
    wavelength2_nm = kwargs.get('wavelength2_nm', 1550)
    dispersion_ps_nm_km = kwargs.get('dispersion_ps_nm_km', 17.0)
    qe1 = kwargs.get('qe1', 0.7)
    dcr1_hz = kwargs.get('dcr1_hz', 100)
    jitter1_s = kwargs.get('jitter1_s', 50e-12)
    dead_time1_s = kwargs.get('dead_time1_s', 50e-9)
    qe2 = kwargs.get('qe2', 0.75)
    dcr2_hz = kwargs.get('dcr2_hz', 120)
    jitter2_s = kwargs.get('jitter2_s', 60e-12)
    dead_time2_s = kwargs.get('dead_time2_s', 45e-9)

    # --- Initial Setup ---
    transmission1 = 10**(-path_loss1_db / 10)
    transmission2 = 10**(-path_loss2_db / 10)
    wavelength_diff_nm = wavelength2_nm - wavelength1_nm
    delay_ps = dispersion_ps_nm_km * distance_km * wavelength_diff_nm
    delay_s = delay_ps * 1e-12

    chunk_duration_s = simulation_time_s / num_chunks
    last_detection_time1 = -np.inf
    last_detection_time2 = -np.inf

    for i in range(num_chunks):
        current_chunk_start_time = i * chunk_duration_s

        # --- Generate events for the current chunk ---
        num_pairs_chunk = int(chunk_duration_s * source_pair_rate_hz)
        pair_intervals = np.random.exponential(1.0 / source_pair_rate_hz, num_pairs_chunk)
        pair_creation_times = np.cumsum(pair_intervals) + current_chunk_start_time
        pair_creation_times = pair_creation_times[pair_creation_times < (i + 1) * chunk_duration_s]

        total_efficiency1 = transmission1 * qe1
        detection_mask1 = np.random.rand(len(pair_creation_times)) < total_efficiency1
        photon_times1 = pair_creation_times[detection_mask1]

        total_efficiency2 = transmission2 * qe2
        detection_mask2 = np.random.rand(len(pair_creation_times)) < total_efficiency2
        photon_times2 = pair_creation_times[detection_mask2]

        num_dark_counts1 = np.random.poisson(dcr1_hz * chunk_duration_s)
        dark_count_times1 = np.random.uniform(current_chunk_start_time, (i + 1) * chunk_duration_s, num_dark_counts1)

        num_dark_counts2 = np.random.poisson(dcr2_hz * chunk_duration_s)
        dark_count_times2 = np.random.uniform(current_chunk_start_time, (i + 1) * chunk_duration_s, num_dark_counts2)

        all_events1 = np.concatenate((photon_times1, dark_count_times1))
        all_events1 += np.random.normal(0, jitter1_s, len(all_events1))
        all_events1.sort()

        all_events2 = np.concatenate((photon_times2, dark_count_times2))
        all_events2 += delay_s
        all_events2 += np.random.normal(0, jitter2_s, len(all_events2))
        all_events2.sort()

        # --- Apply dead time with state ---
        final_tags1, last_detection_time1 = _apply_dead_time_iterative(all_events1, dead_time1_s, last_detection_time1)
        final_tags2, last_detection_time2 = _apply_dead_time_iterative(all_events2, dead_time2_s, last_detection_time2)

        yield final_tags1, final_tags2


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
