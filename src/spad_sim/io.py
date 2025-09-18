import pandas as pd

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
