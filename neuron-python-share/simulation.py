# Simulation and experiment classes / functions

# TODO: Add information on how to run simulations

# ---------------------------------------------------------
# NEURON setup

from neuron import h
h.load_file('stdrun.hoc')

# ---------------------------------------------------------
# Detect if neuron fires

class SpikeDetector:

    def __init__(self, section, v_init=-65, threshold=0, timeout=5):
        # Cell section to monitor (e.g. cell.axon[0] monitors the cell's axon)
        self.section = section 
        
        # Simulation parameters
        self.v_init = v_init # Resting potential (mV)
        self.threshold = threshold # Spike threshold (mV)
        self.timeout = timeout # Simulation duration (ms)

        # Spike information
        self.fired = False
        self.spike_time = None

        # Spike detection NetCon
        self.detector = h.NetCon(section(0.5)._ref_v, None, sec=section)
        self.detector.threshold = threshold
        self.detector.record(self._handle)

    def _handle(self):
        self.fired = True
        self.spike_time = float(h.t)
        h.stoprun = 1

    def reset(self):
        self.fired = False
        self.spike_time = None

    def run(self):
        # Reset simulation
        self.reset()

        # Run simulation
        h.finitialize(self.v_init)
        h.continuerun(self.timeout)

        return self.fired, self.spike_time

# ---------------------------------------------------------
# Determine minimum pulse scaling factor which causes spike

class ThresholdFinder:

    def __init__(self, cell, pulse, detector, max_attempts=100, tol=1e-2, scale_init=-200):

        self.pulse = pulse # RectangularPulse instance; stimulus affecting cell
        self.detector = detector # SpikeDetector instance; to run threshold search

        self.max_attempts = max_attempts # Number of attempts before quitting if threshold not found
        self.tol = tol # Tolerance to stop binary search
        self.scale_init = scale_init # Initial test

        cell.attach_stim(pulse)

    def run(self, display=False):

        attempts = 0

        fired = False
        fire_time = None

        pulse = self.pulse # Need to specify and build new stimuli
        detector = self.detector # Need to detect spikes

        # PART 1: determine (broad) bounds for firing / no firing

        scale_test = self.scale_init # Initial guess

        pulse.set_stim(scale=scale_test) # Update pulse
        
        fired, fire_time = detector.run() # Spike detector simulation
        if display: print(f'scale = {scale_test} | fired = {fired} | fire_time = {fire_time}')

        if fired: # Upper scale bound found

            scale_fire = scale_test

            while True: # Search for lower scale bound

                attempts += 1
                if attempts == self.max_attempts:
                    print('[ThresholdFinder: run] Max number of attempts reached. Search failed.')
                    return

                # Re-run spike detection with new scaling factor
                scale_test = scale_fire / 2
                pulse.set_stim(scale=scale_test)
                fired, fire_time = detector.run()
                if display: print(f'scale = {scale_test} | fired = {fired} | fire_time = {fire_time}')
                
                if fired: 
                    scale_fire = scale_test # Update upper scale bound, try again
                else: 
                    scale_nofire = scale_test # Lower scale bound found
                    break

        else: # Lower scale bound found

            scale_nofire = scale_test

            while True: # Search for upper scale bound

                attempts += 1
                if attempts == self.max_attempts:
                    print('[ThresholdFinder: run] Max number of attempts reached. Search failed.')
                    return

                # Re-run spike detection with new scaling factor
                scale_test = scale_nofire * 2
                pulse.set_stim(scale=scale_test)
                fired, fire_time = detector.run()
                if display: print(f'scale = {scale_test} | fired = {fired} | fire_time = {fire_time}')

                if fired:
                    scale_fire = scale_test # Upper scale bound found
                    break
                else:
                    scale_nofire = scale_test # Update lower scale bound, try again

        # PART 2: binary search for threshold within tolerance

        while abs(scale_fire - scale_nofire) > self.tol:

            attempts += 1
            if attempts == self.max_attempts:
                print('[ThresholdFinder: run] Max number of attempts reached. Search failed.')
                return

            # Run spike detection
            scale_test = (scale_fire + scale_nofire)/2 # Test midpoint
            pulse.set_stim(scale=scale_test)
            fired, fire_time = detector.run()
            if fired: last_fire_time = fire_time
            if display: print(f'scale = {scale_test} | fired = {fired} | fire_time = {fire_time}')

            if fired: scale_fire = scale_test # Decrease upper bound
            else: scale_nofire = scale_test # Increase lower bound

        # Results: Firing threshold, firing time
        self.scale_nofire = scale_nofire
        self.scale_fire = scale_fire
        self.threshold = scale_fire
        self.fire_time = last_fire_time

        return scale_fire
