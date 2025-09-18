# SPAD Simulation

This project is a simulation tool for Single-Photon Avalanche Diode (SPAD) systems, based on the concepts described in https://arxiv.org/abs/2301.08475.

It simulates the time tags of photon arrivals at two detectors, taking into account various physical effects such as quantum efficiency, dark counts, timing jitter, dead time, and chromatic dispersion.

## Installation

1.  Clone the repository:
    ```bash
    git clone <repository-url>
    cd spad-sim
    ```

2.  Install the project and its dependencies using `pip`:
    ```bash
    pip install -e .
    ```

## Usage

The simulation is run from the command line using the `spad-sim` script.

```bash
spad-sim --help
```

This will show all the available parameters you can use to configure the simulation.

### Example

To run a simulation with a 10km fiber, and save the output to `my_simulation.csv`, you can run:

```bash
spad-sim --distance_km 10 --output_file my_simulation.csv
```

The script will print the simulation progress and save the time tags to the specified CSV file. A plot of the coincidence histogram will be displayed by default. You can disable this with the `--no-plot` flag.
