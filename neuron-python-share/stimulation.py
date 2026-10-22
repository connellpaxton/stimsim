# Stimulation classes

# TODO: Secondary electric field effects

# ---------------------------------------------------------
# NEURON setup

from neuron import h

# ---------------------------------------------------------
# EField

class EField:

    def __init__(self, E=[0, 0, 1]):
        self.set_field(E=E)

    def __str__(self):
        return f'[EField] Ex = {self.Ex}, Ey = {self.Ey}, Ez = {self.Ez}, |E| = {self.magnitude}'

    def set_field(self, E):
        self.E = E
        self.Ex = E[0]
        self.Ey = E[1]
        self.Ez = E[2]
        self.magnitude = (E[0]**2 + E[1]**2 + E[2]**2)**0.5

    def normalize(self):
        if self.magnitude == 0: return
        E_normalized = [
            self.Ex / self.magnitude,
            self.Ey / self.magnitude,
            self.Ez / self.magnitude
        ]
        self.set_field(E=E_normalized)

# ---------------------------------------------------------
# Stimulus waveform

class RectangularPulse:

    def __init__(self, delay=0.0, duration=0.1, scale=1):
        # Pulse parameters
        self.delay = delay
        self.duration = duration
        self.scale = scale

        # NEURON Vectors (initialize once)
        self.stim_scale = h.Vector(6)
        self.stim_time = h.Vector(6)

        # Update vectors
        self.set_stim(delay, duration, scale)

    def __str__(self):
        return f'[RectangularPulse] delay = {self.delay}, duration = {self.duration}, scale = {self.scale}'

    def set_stim(self, delay=None, duration=None, scale=None):
        # Update parameters only if provided; else use previous values / defaults
        delay = self.delay if delay is None else delay
        duration = self.duration if duration is None else duration
        scale = self.scale if scale is None else scale

        # Build NEURON Vectors to be attached to a Cell

        # Stimulation vector
        self.stim_scale.x[0] = 0
        self.stim_scale.x[1] = 0
        self.stim_scale.x[2] = scale
        self.stim_scale.x[3] = scale
        self.stim_scale.x[4] = 0
        self.stim_scale.x[5] = 0

        # Time vector
        self.stim_time.x[0] = 0
        self.stim_time.x[1] = delay
        self.stim_time.x[2] = delay
        self.stim_time.x[3] = delay + duration
        self.stim_time.x[4] = delay + duration
        self.stim_time.x[5] = delay + duration + 1
