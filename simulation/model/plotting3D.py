# This has been updated and polished by AI

import numpy as np
import matplotlib.pyplot as plt

from neuron import h
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection


class MembranePotential3D:
    """
    Record and animate NEURON membrane potentials on a 3D morphology.

    Geometry is taken from:
        seg.x_xtra
        seg.y_xtra
        seg.z_xtra

    Voltage is taken from:
        seg.v

    Typical usage
    --------------
    animator = MembranePotential3D(cell.all)

    # Record every 0.1 ms instead of every simulation time step.
    animator.record(record_dt=0.1)

    # Run NEURON simulation here...
    h.finitialize(-70)
    h.continuerun(100)

    animator.load_recording()

    animator.show(
        fps=30,
        playback_speed=0.1,
    )
    """

    def __init__(self, sections):
        self.sections = list(sections)

        if not self.sections:
            raise ValueError("No sections were provided.")

        # ------------------------------------------------------------------
        # Build segment and geometry information once.
        # ------------------------------------------------------------------

        self.segments = [
            seg
            for sec in self.sections
            for seg in sec
        ]

        if not self.segments:
            raise ValueError("No computational segments were found.")

        self._build_geometry()

        # ------------------------------------------------------------------
        # Recording data
        # ------------------------------------------------------------------

        self.t_vec = None
        self.v_vecs = []

        self.t = None
        self.V = None

        # ------------------------------------------------------------------
        # Matplotlib objects
        # ------------------------------------------------------------------

        self.fig = None
        self.ax = None
        self.animation = None
        self.line_collection = None

    # ======================================================================
    # Geometry
    # ======================================================================

    @staticmethod
    def _position(seg):
        """Return the 3D position of a segment."""

        return np.array(
            [
                seg.x_xtra,
                seg.y_xtra,
                seg.z_xtra,
            ],
            dtype=float,
        )

    @staticmethod
    def _segment_key(seg):
        """
        Identify a segment by its section and normalized position.

        Using the section and seg.x avoids depending on the identity of
        the Python wrapper returned for a NEURON segment.
        """

        return seg.sec, float(seg.x)

    def _build_geometry(self):
        """
        Build the static morphology used by the animation.

        Each displayed line is associated with two NEURON segments.
        The voltage of the line is obtained by averaging the voltages
        at those two segments.
        """

        segment_index = {
            self._segment_key(seg): i
            for i, seg in enumerate(self.segments)
        }

        lines = []
        segment_pairs = []

        # --------------------------------------------------------------
        # Lines within each section
        # --------------------------------------------------------------

        section_segments = {
            sec: list(sec)
            for sec in self.sections
        }

        for sec, segments in section_segments.items():

            for seg1, seg2 in zip(segments[:-1], segments[1:]):

                lines.append([
                    self._position(seg1),
                    self._position(seg2),
                ])

                segment_pairs.append([
                    segment_index[self._segment_key(seg1)],
                    segment_index[self._segment_key(seg2)],
                ])

        # --------------------------------------------------------------
        # Lines connecting parent and child sections
        # --------------------------------------------------------------

        included_sections = set(self.sections)

        for child_sec in self.sections:

            ref = h.SectionRef(sec=child_sec)

            if not ref.has_parent():
                continue

            parent_sec = ref.parent

            if parent_sec not in included_sections:
                continue

            parent_segments = section_segments[parent_sec]
            child_segments = section_segments[child_sec]

            if not parent_segments or not child_segments:
                continue

            # Find the parent segment nearest the child connection.
            parent_x = child_sec.parentseg().x

            parent_seg = min(
                parent_segments,
                key=lambda seg: abs(seg.x - parent_x),
            )

            # Find the child segment nearest the connected end.
            if child_sec.orientation() == 0:
                child_seg = child_segments[0]
            else:
                child_seg = child_segments[-1]

            lines.append([
                self._position(parent_seg),
                self._position(child_seg),
            ])

            segment_pairs.append([
                segment_index[self._segment_key(parent_seg)],
                segment_index[self._segment_key(child_seg)],
            ])

        if not lines:
            raise ValueError("No drawable morphology lines were found.")

        self.lines = np.asarray(lines, dtype=float)
        self.segment_pairs = np.asarray(segment_pairs, dtype=int)

    # ======================================================================
    # Recording
    # ======================================================================

    def record(self, record_dt=None, time_step=None):
        """
        Set up NEURON Vector recording.

        Parameters
        ----------
        record_dt : float, optional
            Recording interval in ms.

            For example:
                record_dt=0.1

            records every 0.1 ms.

        time_step : int, optional
            Alternative to record_dt.

            For example:
                time_step=5

            records every 5 * h.dt.

            This is retained for compatibility with the original
            implementation.

        Notes
        -----
        Call this BEFORE h.finitialize() and before running the simulation.
        """

        if record_dt is not None and time_step is not None:
            raise ValueError(
                "Specify either record_dt or time_step, not both."
            )

        if record_dt is None:
            if time_step is None:
                record_dt = h.dt
            else:
                record_dt = time_step * h.dt

        record_dt = float(record_dt)

        if record_dt <= 0:
            raise ValueError("record_dt must be greater than zero.")

        # Time vector
        self.t_vec = h.Vector()
        self.t_vec.record(h._ref_t, record_dt)

        # Voltage vectors
        self.v_vecs = []

        for seg in self.segments:
            vec = h.Vector()
            vec.record(seg._ref_v, record_dt)
            self.v_vecs.append(vec)

        # Clear old NumPy data
        self.t = None
        self.V = None

        self.record_dt = record_dt

    def load_recording(self):
        """
        Convert NEURON Vector recordings into NumPy arrays.

        After this:

            self.t.shape == (n_times,)

            self.V.shape == (n_segments, n_times)
        """

        if self.t_vec is None or not self.v_vecs:
            raise RuntimeError(
                "Recording has not been configured. "
                "Call record() before running the simulation."
            )

        self.t = np.asarray(self.t_vec, dtype=float).copy()

        self.V = np.asarray(
            [np.asarray(vec, dtype=np.float32) for vec in self.v_vecs]
        )

        if len(self.t) == 0:
            raise RuntimeError("No simulation data were recorded.")

        if self.V.shape[1] != len(self.t):
            raise RuntimeError(
                "Voltage and time recordings have different lengths."
            )

        print(f'Sizes: t {self.t.shape}, V {self.V.shape}')

        return self.t, self.V

    # ======================================================================
    # Voltage interpolation
    # ======================================================================

    def _interpolated_voltage(self, time):
        """
        Return interpolated voltage at an arbitrary simulation time.

        The interpolation is performed for every NEURON segment.

        This avoids calling np.interp separately for every segment.
        """

        if time <= self.t[0]:
            return self.V[:, 0]

        if time >= self.t[-1]:
            return self.V[:, -1]

        # Find the recorded samples surrounding this time.
        right = np.searchsorted(self.t, time)

        left = right - 1

        t0 = self.t[left]
        t1 = self.t[right]

        # Fraction between the two recorded samples.
        alpha = (time - t0) / (t1 - t0)

        # Linear interpolation for ALL segments simultaneously.
        return (
            self.V[:, left] * (1.0 - alpha)
            + self.V[:, right] * alpha
        )

    def _line_voltage(self, segment_voltage):
        """
        Convert segment voltages into voltages for displayed morphology lines.

        Each line gets the mean voltage of its two endpoint segments.
        """

        return segment_voltage[self.segment_pairs].mean(axis=1)

    # ======================================================================
    # Animation
    # ======================================================================

    def animate(
        self,
        vlim=(-80, 40),
        cmap="turbo",
        linewidth=2,
        fps=30,
        playback_speed=0.001,
        figsize=(9, 8),
        repeat=True,
    ):
        """
        Create a smooth animation from the recorded data.

        Parameters
        ----------
        vlim : tuple
            Fixed color limits in mV.

        cmap : str
            Matplotlib colormap.

        linewidth : float
            Width of morphology lines.

        fps : float
            Number of animation frames displayed per second.

        playback_speed : float
            Ratio of simulation time to real time.

            1.0  = real time
            0.1  = 10x slower
            0.5  = 2x slower
            2.0  = 2x faster

        figsize : tuple
            Figure size.

        repeat : bool
            Whether to repeat the animation.
        """

        if self.t is None or self.V is None:
            raise RuntimeError(
                "No recording is loaded. "
                "Call load_recording() after the simulation."
            )

        if fps <= 0:
            raise ValueError("fps must be greater than zero.")

        if playback_speed <= 0:
            raise ValueError("playback_speed must be greater than zero.")

        # ------------------------------------------------------------------
        # Figure
        # ------------------------------------------------------------------

        self.fig = plt.figure(figsize=figsize)

        self.ax = self.fig.add_subplot(
            111,
            projection="3d",
        )

        norm = Normalize(
            vmin=vlim[0],
            vmax=vlim[1],
        )

        self.line_collection = Line3DCollection(
            self.lines,
            cmap=cmap,
            norm=norm,
            linewidth=linewidth,
        )

        self.ax.add_collection3d(self.line_collection)

        # Temporary synapse locations TEST
        # syn_x = [795.744, 695.488, 796.032, 695.816]
        # syn_y = [723.648, 768.320, 723.688, 768.336]
        # syn_z = [936.440, 873.640, 936.280, 873.840]
    
        # self.ax.scatter(
        #     syn_x,
        #     syn_y,
        #     syn_z,
        #     color="black",
        #     s=30,
        #     depthshade=False
        # )

        # ------------------------------------------------------------------
        # Set equal spatial limits
        # ------------------------------------------------------------------

        xyz = self.lines.reshape(-1, 3)

        xyz_min = xyz.min(axis=0)
        xyz_max = xyz.max(axis=0)

        center = (xyz_min + xyz_max) / 2
        span = np.max(xyz_max - xyz_min)

        if span <= 0:
            span = 1.0

        half_span = 0.525 * span

        self.ax.set_xlim(
            center[0] - half_span,
            center[0] + half_span,
        )

        self.ax.set_ylim(
            center[1] - half_span,
            center[1] + half_span,
        )

        self.ax.set_zlim(
            center[2] - half_span,
            center[2] + half_span,
        )

        self.ax.set_box_aspect((1, 1, 1))

        self.ax.set_xlabel("x (µm)")
        self.ax.set_ylabel("y (µm)")
        self.ax.set_zlabel("z (µm)")

        # ------------------------------------------------------------------
        # Colorbar
        # ------------------------------------------------------------------

        colorbar = self.fig.colorbar(
            self.line_collection,
            ax=self.ax,
            pad=0.10,
            shrink=0.75,
        )

        colorbar.set_label("Membrane potential (mV)")

        # ------------------------------------------------------------------
        # Animation timing
        # ------------------------------------------------------------------

        # Real-world time between displayed frames, in milliseconds.
        display_interval = 1000.0 / fps

        # Amount of simulation time represented by each animation frame.
        simulation_interval = (
            display_interval * playback_speed
        )

        animation_times = np.arange(
            self.t[0],
            self.t[-1],
            simulation_interval,
        )

        # Make sure the final frame is included.
        if len(animation_times) == 0 or animation_times[-1] < self.t[-1]:
            animation_times = np.append(
                animation_times,
                self.t[-1],
            )

        # ------------------------------------------------------------------
        # Update function
        # ------------------------------------------------------------------

        def update(time):

            # Interpolate voltage between recorded samples.
            segment_voltage = self._interpolated_voltage(time)

            # Convert segment voltages into line voltages.
            line_voltage = self._line_voltage(segment_voltage)

            self.line_collection.set_array(line_voltage)

            self.ax.set_title(
                f"Membrane Potential    t = {time:.2f} ms"
            )

            return self.line_collection,

        # Initialize the first frame.
        update(animation_times[0])

        self.animation = FuncAnimation(
            self.fig,
            update,
            frames=animation_times,
            interval=display_interval,
            blit=False,
            repeat=repeat,
        )

        self.fig.tight_layout()

        return self.fig, self.ax, self.animation

    # ======================================================================
    # Convenience method
    # ======================================================================

    def show(self, **kwargs):
        """Create and display the animation."""

        result = self.animate(**kwargs)

        plt.show()

        return result

    # ======================================================================
    # Geometry update
    # ======================================================================

    def rebuild_geometry(self):
        """
        Rebuild geometry if x_xtra/y_xtra/z_xtra have changed.
        """

        self._build_geometry()
