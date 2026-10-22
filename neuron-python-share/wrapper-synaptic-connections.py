'''
Date  : Last updated 10/22/2026 by Mark
Info  : This script runs a simulation with IClamp stimulation and true synapse connections between N cells.
Input : Specify (1) a list of paths to SWC files, and (2) a path to a synaptic connectivity table CSV file.
Plots : Displays (1) 3D membrane potential animation, and (2) membrane potential vs. time

Notes 
    L*_all.csv files contain ALL synaptic connections in each layer. The SynapseConstructor automatically 
    inserts all synapses from the provided table into the corresponding cells present in the NEURON simulation.

'''

# =========================================================
# USER-DEFINED PARAMETERS

# SWC files for cells
files = [
    'W:/synapse_data_nm/SWC/SWC_downloads_version_1300/SWC_L5/864691135640672571.swc',
    'W:/synapse_data_nm/SWC/SWC_downloads_version_1300/SWC_L5/864691135993429185.swc',
    'W:/synapse_data_nm/SWC/SWC_downloads_version_1300/SWC_L5/864691136025333561.swc'
]

# CSV file for synapses (connectivity table)
connectivity = 'L5_all.csv' 

# Load & show animation and other figure? (True / False)
show_animation = True
show_figure = True

# =========================================================
# Imports 

# Standard libraries
import pandas as pd
import matplotlib.pyplot as plt
from neuron import h
h.load_file('stdrun.hoc')

# Custom Python classes
from cell import Cell, SynapseConstructor
from plotting3D import MembranePotential3D

# ---------------------------------------------------------
# Load cells

print('\nCELLS')

cells = []

for swc in files:
    cell = Cell(swc, nseg_reduction_factor=2)
    cells.append(cell)
    print(cell)

# ---------------------------------------------------------
# Insert synapses 

print('\nSYNAPSES')

sc = SynapseConstructor()
sc.insert_synapses(connectivity)

for pre_syn in sc.pre_syns:
    print(f'Pre-synaptic segment: {pre_syn.preseg()} | Post-synaptic segment: {pre_syn.postseg()}')

# ---------------------------------------------------------
# Set up stimulation

stim = h.IClamp(cells[0].soma[0](0.5))
stim.amp = 1   # nA
stim.delay = 0 # ms
stim.dur = 1   # ms

# ---------------------------------------------------------
# Set up animation

if show_animation:
    plotter = MembranePotential3D(sections=h.allsec()) # Can also specify sections=cells[i].axon, etc.
    plotter.record()

# ---------------------------------------------------------
# Set up vectors for other figure

if show_figure:
    # Time vector
    t = h.Vector().record(h._ref_t)

    # Membrane potential vectors
    v_vecs = []
    for cell in cells:
        v = h.Vector().record(cell.axon[0](0.5)._ref_v)
        v_vecs.append(v)

# ---------------------------------------------------------
# Run simulation

h.dt = 0.001 # ms

print('\nSIMULATION')
print(f'dt = {h.dt}')
print('Running simulation...')

h.finitialize(-65) # mV
h.continuerun(25) # ms

print('Simulation done.')

# ---------------------------------------------------------
# Show figures

if show_animation:
    plotter.load_recording()
    plotter.show(
        vlim=(-80, 40),
        fps=30,
        playback_speed=0.001
    )

if show_figure:
    for v, cell in zip(v_vecs, cells): plt.plot(t, v, label=f'{cell}.axon[0](0.5)')
    plt.title(f'Membrane potential (axon) vs. time for {len(cells)} interconnected cells')
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane potential (mV)')
    plt.legend(loc='upper right')
    plt.show()