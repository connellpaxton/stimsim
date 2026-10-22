# Cell constructor classes

'''
TODO: Passable membrane parameters
  - Cell class argument: folder containing mod files with naming conventon for sections (e.g., Ca_axon.mod)
  - If no mech folder provided, use default folder "mech" or add mechanisms as it is now

TODO: Automatic mechanism insertion 
  - In mechanism folder, mod files can be named based on the section(s) they apply to,
    e.g., Ca_soma.mod, K_soma_axon_dend.mod, xtra_all.mod

TODO: Add MultipleCells class for efficiently creating many cells
  - Input could be list of file names, directory name with files in it, or connectivity table w/ file locations specified

TODO: Vectorize extracellular potential computation

TODO: Investigate and fix mech file problems

TODO: Update links in SynapseConstructor.update_params & .__init__() when moving over to NEURON 9.0.2

'''

# ---------------------------------------------------------
# NEURON setup

from neuron import h
h.load_file('import3d.hoc')

from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------
# Cell setup

class Cell:

    def __init__(self, file, nseg_reduction_factor=0):
        '''
        Builds a cell in NEURON from the provided SWC file.

        Inputs:
            file (str): SWC file name.
            nseg_reduction_factor (int >= 0): Specify to reduce nseg, lowering computational load. Recommended range: 0-3.
        
        '''
        # Set class attributes
        self.file = file
        self.morphology = Path(file).stem
        self.nrf = max(int(nseg_reduction_factor), 0)

        # Build NEURON cell
        self._load_morphology() # Build cell from SWC file
        self._add_biophysics() # Add membrane mechanisms
        self._discretize() # Set discretization (nseg)
        self._interpolate_nodes() # Get computational nodes
        self._set_pointers() # Set up xtra mechanism

    def __str__(self):
        # What to print when referring to object
        return f'{self.morphology}'

    def _load_morphology(self):
        # Load cell morphology from SWC file
        cell = h.Import3d_SWC_read()
        cell.input(f'{self.file}')
        i3d = h.Import3d_GUI(cell, False)
        i3d.instantiate(self)

    def _add_HH_biophysics(self): # Temporary; for testing purposes
        # Add channel parameters
        for sec in self.all:
            sec.Ra = 100
            sec.cm = 1.0
            sec.insert('hh')
            sec.insert('extracellular')
            sec.insert('xtra')

    def _add_biophysics(self):
        # Add channel parameters (following Neuron_Resources/864691134918921738/parameters.hoc)

        if hasattr(self, 'all'):
            for sec in self.all:
                sec.Ra = 100
                sec.cm = 1
                sec.insert('pas')
                for seg in sec:
                    seg.pas.g = 6.3e-5
                    seg.pas.e = -64.601696
                sec.insert('extracellular')
                sec.insert('xtra') 

        if hasattr(self, 'axon'):
            for sec in self.axon:
                sec.insert('Ca')
                sec.insert('CaDynamics_E2')
                sec.insert('Ca_LVAst')
                sec.insert('Im')
                sec.insert('K_Pst') # effect: repeated firing
                sec.insert('K_Tst')
                sec.insert('NaTa_t') # effect: no spike, slowly decreasing potential
                sec.insert('Nap_Et2')
                sec.insert('SK_E2')
                sec.insert('SKv3_1') # effect: weird oscillations, stablizes at -30 mV

        if hasattr(self, 'soma'):
            for sec in self.soma:
                for seg in sec:
                    seg.pas.g = 2e-5
                    seg.pas.e = -69.781406
                sec.insert('Ca_soma')
                sec.insert('CaDynamics_E2_soma')
                sec.insert('Ca_LVAst_soma')
                sec.insert('Im_soma')
                sec.insert('K_Pst_soma')
                sec.insert('K_Tst_soma')
                sec.insert('NaTs2_t_soma') # effects: decreases spike height
                sec.insert('Nap_Et2_soma')
                sec.insert('SK_E2_soma')
                sec.insert('SKv3_1_soma') # effects: increases lingering, repeated firing

        if hasattr(self, 'dend'):
            for sec in self.dend:
                for seg in sec:
                    seg.pas.g = 1e-6 # Super low?
                    seg.pas.e = -63.118492
                sec.insert('Ih')
                sec.insert('Im_apic')
                sec.insert('K_Pst_apic')
                sec.insert('K_Tst_apic')
                sec.insert('NaTs2_t_apic')
                sec.insert('Nap_Et2')
                sec.insert('SKv3_1_apic')

        if hasattr(self, 'apic'):
            for sec in self.apic:
                for seg in sec:
                    seg.pas.g = 1e-6
                    seg.pas.e = -63.118492
                sec.insert('Ih')
                sec.insert('Im_apic')
                sec.insert('K_Pst_apic')
                sec.insert('K_Tst_apic')
                sec.insert('NaTs2_t_apic')
                sec.insert('Nap_Et2')
                sec.insert('SKv3_1_apic')

        if hasattr(self, 'dend_5'):
            for sec in self.dend_5:
                for seg in sec:
                    seg.pas.g = (8/9)*1e-6
                sec.cm = 0.02

        if hasattr(self, 'dend_6'):
            for sec in self.dend_6:
                sec.insert('K_Pst')
                sec.insert('K_Tst')
                sec.insert('NaTa_t_nor')
                sec.insert('Nap_Et2')
                sec.insert('SKv3_1')

    def _add_biophysics_alternate(self): # Temporary; for testing purposes
        # Add channel parameters (following 00__Neuron_Example/parameters.hoc)
        
        if hasattr(self, 'all'):
            for sec in self.all:
                sec.Ra = 100
                sec.insert('pas')
                for seg in sec:
                    seg.pas.g = 3e-5
                    seg.pas.e = -75
                sec.insert('extracellular')
                sec.insert('xtra')

        # Mechanisms for axon sections
        if hasattr(self, 'axon'):
            for sec in self.axon:
                sec.cm = 1
                sec.insert('CaDynamics_E2')
                sec.insert('Ca_HVA')
                sec.insert('Ca_LVAst')
                sec.insert('K_Pst')
                sec.insert('K_Tst')
                sec.insert('NaTa_t')
                sec.insert('Nap_Et2')
                sec.insert('SK_E2')
                sec.insert('SKv3_1')

        # Mechanisms for soma sections
        if hasattr(self, 'soma'):
            for sec in self.soma:
                sec.cm = 1
                sec.insert('CaDynamics_E2_soma')
                sec.insert('Ca_HVA_soma')
                sec.insert('Ca_LVAst_soma')
                sec.insert('Ih')
                sec.insert('NaTs2_t_soma')
                sec.insert('SK_E2_soma')
                sec.insert('SKv3_1_soma')

        # Mechanisms for dendrite sections
        if hasattr(self, 'dend'):
            for sec in self.dend:
                sec.cm = 2
                sec.insert('Ih')

        # Mechanisms for apical sections
        if hasattr(self, 'apic'):
            for sec in self.apic:
                sec.cm = 2
                sec.insert('Ih')
                sec.insert('Im')
                sec.insert('NaTs2_t_apic')
                sec.insert('SKv3_1_apic')

    def _discretize(self):
        '''
        Sets discritization (nseg) using "d_lambda" method. Needs to happen AFTER adding mechanisms
        because h.d_lambda() depends on inserted mechanisms.

        Reduction of nseg up to a factor of ~3 yields sub-percent agreement with max nseg results.
        
        '''
        if hasattr(self, 'all'):
            freq = 100 # Frequency (Hz) at which AC length constant will be computed | tutorial default = 100
            d_lambda = 0.01
            for sec in self.all:
                # sec.push() and h.pop_section() ensure current section is being accessed
                sec.push()
                nseg = int((sec.L/(d_lambda * h.lambda_f(freq)) + 0.999)/2)*2 + 1 # "d_lambda" method
                nseg = int(nseg / (3**self.nrf) / 2)*2 + 1 # nseg reduction
                sec.nseg = nseg
                h.pop_section()

    def _interpolate_nodes(self):
        # Compute x,y,z coordinates of computational nodes in model cell from pt3d data
        # (Formerly "grindaway()" from "interpxyz.hoc"; code and comments taken from there)
        
        if hasattr(self, 'all'):
            for sec in self.all:
                if h.ismembrane('xtra', sec=sec):
                    # Get data for the section
                    nn = sec.n3d()
                    xx = h.Vector(nn)
                    yy = h.Vector(nn)
                    zz = h.Vector(nn)
                    length = h.Vector(nn)

                    for ii in range(nn):
                        xx.x[ii] = sec.x3d(ii)
                        yy.x[ii] = sec.y3d(ii)
                        zz.x[ii] = sec.z3d(ii)
                        length.x[ii] = sec.arc3d(ii)

                    # To use Vector class's interpolate() method, must first scale independent variable
                    # e.g., normalize length along centroid 
                    length.div(length.x[nn-1])

                    # Initialize destination "independent vector"
                    rnge = h.Vector(sec.nseg+2)
                    rnge.indgen(1/sec.nseg)
                    rnge.sub(1/(2*sec.nseg))
                    rnge.x[0] = 0
                    rnge.x[sec.nseg+1] = 1

                    # "length" contains the normalized distances of the pt3d points along the centroid
                    #   of the section. These are spaced at irregular intervals.
                    # "rnge" contains the normalized distances of the nodes along the centroid of the 
                    #   section. These are spaced at regular intervals. Ready to interpolate.

                    xint = h.Vector(sec.nseg+2)
                    yint = h.Vector(sec.nseg+2)
                    zint = h.Vector(sec.nseg+2)
                    xint.interpolate(rnge, length, xx)
                    yint.interpolate(rnge, length, yy)
                    zint.interpolate(rnge, length, zz)

                    # For each node, assign x,y,z values to x_xtra, y_xtra, z_xtra
                    # Don't bother computing coordinates of the 0 and 1 ends 
                    # Also avoid writing coorinates of the 1 end into the last interval node's coordinates
                    for ii in range(1, sec.nseg+1):
                        xr = rnge.x[ii]
                        sec(xr).xtra.x = xint.x[ii]
                        sec(xr).xtra.y = yint.x[ii]
                        sec(xr).xtra.z = zint.x[ii]

    def _set_pointers(self):
        # Set pointers (so 'xtra' works) after sections have been created
        if hasattr(self, 'all'):
            for sec in self.all:
                if (h.ismembrane('xtra', sec=sec) and h.ismembrane('extracellular', sec=sec)):
                    for seg in sec:
                        h.setpointer(seg._ref_e_extracellular, 'ex', seg.xtra)

    def compute_extracellular_potential(self, E):
        # E = [Ex, Ey, Ez]
        Ex = E[0]
        Ey = E[1]
        Ez = E[2]

        # Compute extracellular potential for all sections
        if hasattr(self, 'all'):
            for sec in self.all:
                for seg in sec:
                    # Note that xtra positions are in um -> convert to m for calculation
                    ve = -((seg.xtra.x * 1e-6 * Ex) + (seg.xtra.y * 1e-6 * Ey) + (seg.xtra.z * 1e-6 * Ez)) * 1e3
                    seg.xtra.es = ve

    def attach_stim(self, pulse):
        # pulse: RectangularPulse object with attributes:
        #   stim_scale: NEURON Vector representing pulse amplitude
        #   stim_time: NEURON Vector representing pulse on/off time

        stim_scale = pulse.stim_scale
        stim_time = pulse.stim_time

        attached = 0
        if hasattr(self, 'all'):
            for sec in self.all:
                if attached == 0: # Only run if it hasn't been attached yet
                    if h.ismembrane('xtra', sec=sec):
                        stim_scale.play(h._ref_stim_xtra, stim_time, 1)
                        attached = 1

# ---------------------------------------------------------
# Synapse constructor

class SynapseConstructor:

    def __init__(self, post_tau1=0.5, post_tau2=3.0, post_e=0.0, pre_threshold=0.0, pre_delay=1.0, pre_weight=0.01):
        '''
        Initialize synaptic parameters.
        
        Inputs:
            post_tau1 (float): Rise time constant; when synaptic event arrives, rises to max conductance 
            post_tau2 (float): Decay time constant; controls duration of synaptic event
            post_e (float): Reversal potential (direction of potential change); controls excitatory/inhibitory behavior
            pre_threshold (float): Membrane potential at which synaptic event is triggered.
            pre_delay (float): Delay (ms) after which to send event signal to target (e.g., an Exp2Syn); simulates prop. delay
            pre_weight (float): Strength of synaptic event; weight = peak conductance (uS); controls synaptic current
        
        Note: See the following links for more information.
            NetCon: https://www.neuronsimulator.org/en/8.2.3/python/modelspec/programmatic/network/netcon.html#netcon
            Exp2Syn: https://www.neuronsimulator.org/en/8.2.3/python/modelspec/programmatic/mechanisms/mech.html#Exp2Syn

        '''
        self.update_params(
            post_tau1=post_tau1, 
            post_tau2=post_tau2,
            post_e=post_e, 
            pre_threshold=pre_threshold, 
            pre_delay=pre_delay, 
            pre_weight=pre_weight
        )

    @staticmethod
    def find_synapse_segment(cell, synapse_xyz, units='nm'):
        '''
        Determines and returns the segment closest in space to the provided synapse coordinates.
    
        Inputs:
            cell (Cell object): NEURON cell which is known to have the synapse described by synapse_xyz.
            synapse_xyz (list or numpy array): [x, y, z] coordinates of synapse on cell.
            units (str): 'nm' (synapse table default) or 'um' (NEURON & swc default); scales input if units='nm'
    
        Outputs:
            synapse_seg (NEURON segment): Segment belonging to cell which is most closely located to the synapse.
        
        '''
        # Build lists containing segment information
    
        segment_list = []
        coordinates_list = []
    
        for sec in cell.all:
            for seg in sec:
                segment_list.append(seg)
                x = seg.x_xtra
                y = seg.y_xtra
                z = seg.z_xtra
                coordinates_list.append([x, y, z])
    
        # Convert numerical lists to numpy for efficient search
    
        coordinates = np.array(coordinates_list) # Matrix of size (N_segments, 3)
        synapse_xyz = np.array(synapse_xyz) # Vector of length 3
        if units=='nm': synapse_xyz = synapse_xyz / 1e3 # Convert to um if originally nm
    
        # Find index where coordinates most closely matches synapse_xyz
    
        distances = np.linalg.norm(coordinates - synapse_xyz, axis=1)
        idx = distances.argmin()
        synapse_seg = segment_list[idx]
    
        return synapse_seg

    @staticmethod
    def get_all_cells():
        '''
        Returns two lists: one containing all cell objects, and one containing all cell names.
    
        Inputs:
            None. Implements h.allsec() to identify cells from their sections.
    
        Outputs:
            cell_objects (list of objects): list of cell objects, e.g., [cell1, cell2]
            cell_names (list of strings): list of cell names, e.g., ['864691134886404090', '864691134886653434']
        
        '''
        # Initialize set to store unique cell objects
        cell_objects = set()
    
        # For each section, determine the parent cell (NEURON has no other way to get all cells)
        for sec in h.allsec():
            cell = sec.cell()
            if cell is not None:
                cell_objects.add(cell)
    
        # Turn into list; make list of cell names
        cell_objects = list(cell_objects)
        cell_names = [str(c) for c in cell_objects]
    
        return cell_objects, cell_names 

    def update_params(self, post_tau1=None, post_tau2=None, post_e=None, pre_threshold=None, pre_delay=None, pre_weight=None):
        '''
        Modify synaptic parameters; run before insert_synapses.
        
        Inputs:
            post_tau1 (float): Rise time constant; when synaptic event arrives, rises to max conductance 
            post_tau2 (float): Decay time constant; controls duration of synaptic event
            post_e (float): Reversal potential (direction of potential change); controls excitatory/inhibitory behavior
            pre_threshold (float): Membrane potential at which synaptic event is triggered.
            pre_delay (float): Delay (ms) after which to send event signal to target (e.g., an Exp2Syn); simulates prop. delay
            pre_weight (float): Strength of synaptic event; weight = peak conductance (uS); controls synaptic current
        
        Note: See the following links for more information.
            NetCon: https://www.neuronsimulator.org/en/8.2.3/python/modelspec/programmatic/network/netcon.html#netcon
            Exp2Syn: https://www.neuronsimulator.org/en/8.2.3/python/modelspec/programmatic/mechanisms/mech.html#Exp2Syn
        
        '''

        # Check: post_tau2 > post_tau1? (Required by NEURON)
        if post_tau2 <= post_tau1: raise ValueError('post_tau2 must be greater than post_tau1.')
        
        # Post-synaptic parameters (only update if given as function argument; else keep old value)
        self.tau1 = self.tau1 if post_tau1 is None else post_tau1
        self.tau2 = self.tau2 if post_tau2 is None else post_tau2
        self.e = self.e if post_e is None else post_e

        # Pre-synaptic parameters (only update if given as function argument; else keep old value)
        self.threshold = self.threshold if pre_threshold is None else pre_threshold
        self.delay = self.delay if pre_delay is None else pre_delay
        self.weight = self.weight if pre_weight is None else pre_weight

    def insert_synapses(self, connectivity_table):
        '''
        Inserts all synapses described in connectivity_table to the cells which exist in NEURON.

        Inputs:
            connectivity_table (string): csv file containing synaptic connections.
            TODO: cells (list) or ('all'): list of cells in which to insert synapses; default is all cells.

        Outputs: 
            post_syns (list): list of post-synaptic processes (Exp2Syn objects).
            pre_syns (list): list of pre-synaptic processes (NetCon objects).
        
        '''
        # Load connectivity table
        synapse_table = pd.read_csv(
            connectivity_table,
            dtype={
                'pre_pt_root_id': str,
                'post_pt_root_id': str
            }
        )

        # Get list of all cells currently in NEURON model
        cell_objects, cell_names = self.get_all_cells()

        # Initialize container lists
        post_syns = []
        pre_syns = []

        # Loop through each synaptic connection in table
        for _, row in synapse_table.iterrows():

            # Check if cells from table exist in the NEURON model
            if (str(row['pre_pt_root_id']) in cell_names) and (str(row['post_pt_root_id']) in cell_names):

                # Get the cells referenced in the current table row
                for cell_object, cell_name in zip(cell_objects, cell_names):
                    if str(row['post_pt_root_id']) == cell_name: post_cell = cell_object
                    if str(row['pre_pt_root_id']) == cell_name: pre_cell = cell_object
                
                # Load synapse parameters
                post_xyz = list(row[["post_pt_position_x", "post_pt_position_y", "post_pt_position_z"]])
                pre_xyz = list(row[["pre_pt_position_x", "pre_pt_position_y", "pre_pt_position_z"]])
                post_seg = self.find_synapse_segment(post_cell, post_xyz, units='nm')
                pre_seg = self.find_synapse_segment(pre_cell, pre_xyz, units='nm')

                # Insert synapse mechanisms
                post_syn = h.Exp2Syn(post_seg)
                pre_syn = h.NetCon(pre_seg._ref_v, post_syn, sec=pre_seg.sec)

                # Set post-synaptic parameters
                post_syn.tau1 = self.tau1
                post_syn.tau2 = self.tau2
                post_syn.e = self.e

                # Set pre-synaptic parameters
                pre_syn.threshold = self.threshold
                pre_syn.delay = self.delay
                pre_syn.weight[0] = self.weight

                # Append synapses to container lists
                post_syns.append(post_syn)
                pre_syns.append(pre_syn)

        # Save connections
        self.post_syns = post_syns
        self.pre_syns = pre_syns

        # Warn if no synapses were found
        if len(post_syns) == 0 and len(pre_syns) == 0:
            print('SynapseConstructor warning: no synaptic connections were found.')
