import numpy as np
from spad_sim.core import simulate_spad_time_tags

def test_simulation_runs():
    """
    Test that the simulation runs without errors and returns the correct types.
    """
    tags1, tags2 = simulate_spad_time_tags(simulation_time_s=0.01)
    assert isinstance(tags1, np.ndarray)
    assert isinstance(tags2, np.ndarray)
