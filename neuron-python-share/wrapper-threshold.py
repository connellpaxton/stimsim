# Sets up and runs simulation for detecting firing threshold for one cell

# ---------------------------------------------------------
# User-defined parameters

# Cell morphology (.swc file)
file = 'interp_final_typerad.swc'

# E-field components (V/m)
Ex = 0
Ey = 0
Ez = 1

# ---------------------------------------------------------
# NEURON setup

# NEURON itself
from neuron import h
h.load_file('stdrun.hoc')

# Python support classes
from cell import Cell
from stimulation import EField, RectangularPulse
from simulation import SpikeDetector, ThresholdFinder

# ---------------------------------------------------------
# Miscellaneous imports

import pandas as pd # For saving data
import time # For tracking execution speed

# ---------------------------------------------------------
# Cell and stimulation setup (using Python classes)

# Construct cell
cell = Cell(file, nseg_reduction_factor=2)

# Set E-field
field = EField(E=[Ex, Ey, Ez])
field.normalize()

# Compute extracellular potential from E-field
cell.compute_extracellular_potential(field.E)

# Create stimulation pulse 
pulse = RectangularPulse()

# Configure global simulation parameters
h.dt = 0.001

start_time = time.time() # TEST

# Threshold detector simulation
detector = SpikeDetector(cell.axon[0])
finder = ThresholdFinder(cell, pulse, detector)
threshold = finder.run(display=True)

print(f'\nThreshold: {threshold:.4f}\n') 

execution_time = time.time()-start_time # TEST
print(f'Execution time: {execution_time/60:.2f}\n') # TEST

# ---------------------------------------------------------
# Save results

# output = pd.DataFrame([{
#     'threshold': finder.threshold,
#     'fire_time': finder.fire_time,
#     'threshold_upper_bound': finder.scale_fire,
#     'threshold_lower_bound': finder.scale_nofire,
#     'execution_time': execution_time # TEST
# }])
# output.to_csv('threshold.csv', index=False)